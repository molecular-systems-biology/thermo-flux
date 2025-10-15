#! 

""" Module provides function to load and save ThermoModel objects to and from Yaml files. 
Here, the logic is to make thermomodel objects and dict interconvertible, and then use yaml dump/load functions"""

###import packages
import cobra
import cobra.io as cio
from collections import OrderedDict
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING, Dict, List, Sequence, Set, Union
from thermo_flux.core.model import ThermoModel
from thermo_flux.core.reaction import ThermoReaction
from thermo_flux.core.metabolite import ThermoMetabolite
import yaml
from ruamel.yaml.compat import StringIO
from pathlib import Path
from equilibrator_api import ComponentContribution
from equilibrator_assets.local_compound_cache import LocalCompoundCache

YAML_SPEC = "1.2"

### define thermo objects attributes 
##model attributes 

thermo_model_attributes=["pH", "I", "pMg", "T", "phi", "rmse_inf", "max_drG", "phi_dict"]

##reaction attributes
##make a list of all attributes in reaction class
thermo_reaction_attributes=["_drG0", "_drG0prime", "_drG_SE", "_drGtransport", "_drGtransform", "_transported_h", "_transported_charge", "_transported_mets", "_balanced", "_ignore_snd"]

##metabolite attributes
##make a list of all attributes in metabolite class
thermo_metabolite_attributes=["_dfG0", "_dfG0prime", "_dfG_SE", "_redox", "_biomass", "_ignore_conc"]

 

def tf_update_optional(
    cobra_object: "Object",
    new_dict: Dict,
    ordered_keys: Sequence,
) -> None:
    """modified from cobra.io._update_optional
    Update `new_dict` with optional attributes from `cobra_object`.

    Parameters
    ----------
    cobra_object : cobra.Object
        The cobra Object to update optional attributes from.
    new_dict : dict
        The dictionary to update optional attributes for.
    optional_attribute_dict : dict
        The dictionary to use as default value lookup store.
    ordered_keys : list, tuple
        The keys to update values for.

    Raises
    ------
    IndexError
        If key in `ordered_keys` is not found in `optional_attribute_dict`.
    AttributeError
        If key in `ordered_keys` is not found in `cobra_object`.

    """
    for key in ordered_keys:
        value = getattr(cobra_object, key)
        if value is None :
            continue
        new_dict[key] = cio.dict._fix_type(value)


### saving the 

## first i do save model to dict
def thermo_model_to_dict(model: ThermoModel) -> Dict:
    """modified from cobra.io.model_to_dict"""
    #first dump all cobra attributes in a dictionnary
    obj=cio.model_to_dict(model)
    tf_update_optional(model, obj, thermo_model_attributes)

    ##first add model atttributes
    for rxn_dict_i in obj['reactions']: #iterate over all reactions
        for key in thermo_reaction_attributes: #add all rxn attributes
            rxn_dict_i[key]= getattr(model.reactions.get_by_id(rxn_dict_i['id']),key)

    #then metabolites
    for met_dict_i in obj['metabolites']:
        for key in thermo_metabolite_attributes:
            met_dict_i[key]= getattr(model.metabolites.get_by_id(met_dict_i['id']),key)

    return obj
    
#then replicate save function of cobra to dump in yaml
def save_thermo_model(model: ThermoModel, filename: str,**kwargs) -> None:
    """Save a ThermoModel to a YAML file.

    Parameters
    ----------
    model : ThermoModel
        The model to save.
    filename : str
        The filename to save the model to.

    """
    obj = thermo_model_to_dict(model)
    obj["version"] = YAML_SPEC
    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    with open(filename, "w") as file_handle:
            yaml.dump(obj, file_handle,**kwargs )



### loading the model : functions

def drop_metabolite_keys(met: Dict) -> Dict:
    """Drop keys from a metabolite dict that are not in the cobra metabolite class."""
    return {key: val for key, val in met.items() if key not in thermo_metabolite_attributes}

def drop_reaction_keys(rxn: Dict) -> Dict:
    """Drop keys from a reaction dict that are not in the cobra reaction class."""
    return {key: val for key, val in rxn.items() if key not in thermo_reaction_attributes}

def thermo_metabolite_att_from_dict(cobra_met,met_dict: Dict) -> ThermoModel:
    tmet=ThermoMetabolite(cobra_met)
    for key,val in {k:v for k,v in met_dict.items() if k in thermo_metabolite_attributes}.items():
        setattr(tmet,key,val)




###reactions functions 


def load_cache_tmodel(model):#first define a cache, with a given file path if it exists
    lc = LocalCompoundCache()

    cache_name = str(model.id)+'_compound.sqlite'
    from pathlib import Path
    #try and load a local cache with the model id
    if Path(cache_name).is_file(): 
        lc.load_cache(cache_name)
        model._lc = lc
    else: #otherwise generate a new local cache 
        lc.generate_local_cache_from_default_zenodo(cache_name, force_write=True)
        lc.load_cache(cache_name)
        model._lc=lc


def add_thermo_reaction_from_dict(rxn_dict: Dict,model) -> ThermoModel:
    """here i first make a cobra reaction from the dict, then add it to the model. As the model already has a cc attribute and the reaciton is linked to it
    I can turn the reaction in a thermoReaction. Finally i replace the cobra rxn by the thermoreaction"""
    cobra_dict=drop_reaction_keys(rxn_dict)
    cobra_rxn= cio.dict._reaction_from_dict(cobra_dict,model)
    model.add_reactions([cobra_rxn])
    ###initialise the metabolites compound attribute by activate the compound setter
    for met, stoich in model.reactions.get_by_id(rxn_dict['id']).metabolites.items(): #initialisation will be slightly longer but it avoids to compute it each time net_elements is called
        print(met.compound)
    trxn=ThermoReaction(model.reactions.get_by_id(rxn_dict['id']))
    for key,val in {k:v for k,v in rxn_dict.items() if k in thermo_reaction_attributes}.items(): # now add only tf attributes
        setattr(trxn,key,val)
    model.remove_reactions([cobra_rxn])
    model.add_reactions([trxn])


###whole function 
def thermo_model_from_dict(obj: Dict) -> ThermoModel:
    """modified from cobra.io.model_from_dict"""
    if "reactions" not in obj:
        raise ValueError("Object has no .reactions attribute. Cannot load.")
    model=cobra.core.model.Model()

    ##metabolites : add them to cobra model, then add tf attributes 
    model.add_metabolites(
    [ThermoMetabolite(cio.dict._metabolite_from_dict(metabolite)) for metabolite in obj["metabolites"]]
    )
    for metabolite_dict in obj["metabolites"]:
        cobra_met=model.metabolites.get_by_id(metabolite_dict['id'])
        thermo_metabolite_att_from_dict(cobra_met,metabolite_dict)

    ### then i need to connect a component contribution object to the model, so that we can fetch compounds for metabolites when adding rxns
    load_cache_tmodel(model) ##first load cache
    
    model.cc= ComponentContribution(rmse_inf=obj['rmse_inf'], ccache=model._lc.ccache)

    ##add reactions 
    # ### then i  need to add reactions to the model
    for rxn_dict in obj['reactions']:
        add_thermo_reaction_from_dict(rxn_dict,model)

    ##finally transform into a ThermoModel and add attributes from the dict
    tmodel=ThermoModel(model)
    objective_reactions = [
        rxn for rxn in obj["reactions"] if rxn.get("objective_coefficient", 0) != 0
    ]
    coefficients = {
        model.reactions.get_by_id(rxn["id"]): rxn["objective_coefficient"]
        for rxn in objective_reactions
    }
    cobra.util.solver.set_objective(model, coefficients)
    #finally add model attributes and tf specific aatributes
    for k, v in obj.items():
        if k in {"id", "name", "notes", "compartments", "annotation"} or k in thermo_model_attributes and k !='phi_dict':
            setattr(tmodel, k, v)

    return tmodel


def thermo_model_from_yaml(filename: str | Path) -> ThermoModel:
    """Load a ThermoModel from a YAML file.

    Parameters
    ----------
    filename : str
        The filename to load the model from.

    Returns
    -------
    ThermoModel
        The loaded model.

    """
    if isinstance(filename, (str, Path)):
        with open(filename, "r") as file_handle:
            return thermo_model_from_dict(yaml.load(file_handle))
    else:
        return thermo_model_from_dict(yaml.load(filename))