#! 

# Packages needed:
import cobra
import h5py
import pandas as pd

def generate_report(tmodel):
    '''generate a df with reporting on balancing and transporters. Used in update_thermo_info'''
    report_df=pd.DataFrame(index=[rxn.id for rxn in tmodel.reactions], columns=['stoichiometry','balanced','transporter','transported_mets','transported_h','transported_charge','priority'])
    report_df['stoichiometry']=[rxn.reaction for rxn in tmodel.reactions]
    report_df['balanced']=[rxn.balanced for rxn in tmodel.reactions]
    report_df['transporter']=[True  if len(rxn.compartments) > 1 else False for rxn in tmodel.reactions]
    report_df['transported_mets']=[rxn.transported_mets for rxn in tmodel.reactions]
    report_df['transported_h']=[rxn.transported_h for rxn in tmodel.reactions]
    report_df['transported_charge']=[rxn.transported_charge for rxn in tmodel.reactions]
    report_df['priority']=[hasattr(r,'report_prio') for r in tmodel.reactions]
    return report_df.sort_values(by='priority',ascending=False)

def write_model(model, sbml_filename, hdf5_filename=None, attrs_for_hdf5=None, **kwargs):
    """
    Saves the model to an .sbml file; optionally stores additional
    information in a .hdf5 file (nullspace matrix, stoichiometric matrix, etc.)
    """
    write_to_sbml(model, sbml_filename, **kwargs)

    if hdf5_filename is not None:
        write_to_hdf5(model, filename=hdf5_filename, attrs_to_save=attrs_for_hdf5)

    return


# TO-DO: if dealing with pd.DataFrames, may use DataFrame.to_hdf() instead
def write_to_hdf5(model, filename, attrs_to_save):
    """
    Writes additional model information to a .hdf5 file.
    """
    with h5py.File(filename, "w") as f:
        model_grp = f.create_group(model.id)  # /model

        data_grp = model_grp.create_group("data")  # /model/data
        for k in attrs_to_save:
            attr = getattr(model, k)
            values = attr.magnitude  # for objects of Pint
            units = format(attr.units)

            grp = data_grp.create_group(k)  # /model/data/grp
            grp.create_dataset(k, data=values)
            # grp.attrs.create("index", [x for x in values.index])
            # grp.attrs.create("columns", [x for x in values.columns])
            grp.attrs.create("units", units)

    return


def write_to_sbml(model, filename, **kwargs):
    """
    Writes the instance of cobra.Model (or thermo_flux.ThermoModel ???)
    to an .sbml file.

    Uses the appropriate COBRA function to write most of the model information.

    Additional parameters and variables written as "notes" at either the
    model-, reaction-, or metabolite-levels.
    """
    model = arrange_notes(model)

    cobra.io.write_sbml_model(model, filename, **kwargs)

    return


def arrange_notes(model):
    """
    Assigns dictionaries to the instances of either Model, Reaction or
    Metabolite containing the additional information to be stored.
    """
    model = _arrange_metabolite_notes(model)

    model = _arrange_reaction_notes(model)

    model = _arrange_model_notes(model)

    return model


def _arrange_model_notes(model):
    """
    Sets the dictionary with 'parameters' as the 'notes'
    attribute of a cobra.Model
    """

    attrs_to_save = ["I", "pH", "pMg", "phi", "T",
                     "gdiss_lim", "rmse_inf", "max_drG"]
    update_notes(model, attrs_to_save)

    return model


def _arrange_reaction_notes(model):
    """
    Sets the values of Gibbs energy and related thermodynamic information
     to the 'notes' attribute of a cobra.Reaction
    """
    for r in model.reactions:

        attrs_to_save = ["drG0", "drG0prime", "drGtransport", 
                         "drGtransform", "drG_h_transport", "drG_c_transport",
                         "drG", "drG_SE", "transported_h", "transported_charge",
                         "balanced"]
        update_notes(r, attrs_to_save)

    return model


def _arrange_metabolite_notes(model):
    """
    Sets the values of lncLimits to the 'notes'
    attribute of a cobra.Metabolite
    """
    for m in model.metabolites:

        attrs_to_save = ["concentration", "lower_bound", "upper_bound",
                         "accession", "redox", "biomass",
                         "dfG0", "dfG0prime", "dfG_SE"]
        update_notes(m, attrs_to_save)

    return model


def update_notes(obj, attrs_to_save):
    """
    Stores the attributes listed in attrs_to_save in the 'notes' field of 
    each object (Metabolite, Reaction or Model).
    """
    for attr_name in attrs_to_save:
        attr = getattr(obj, attr_name)
        flattened_dict = _flatten_parameters(attr, attr_name)
        obj.notes.update(flattened_dict)
    return


def _flatten_parameters(values, name):
    """
    Returns a dictionary that will be used to populate the 'notes' field of 
    an object.
    Since these must be single-layered dictionaries, this function flattens 
    dictionaries of dictionaries.
    """
    if values is None:
        return {}
    elif isinstance(values, dict):
        return {f"{name}_{key}": f"{value}" for key, value in values.items()}
    else:
        return {f"{name}": f"{values}"}


# LEGACY function
def _helper_set_notes(model, attr_name, obj):
    """
    Populates the 'notes' attribute of object 'obj' (cobra.Reaction or
    cobra.Metabolite) with the data taken from the appropriate model
    attribute ('attr_name').
    Data is selected on the basis of the object's ID.
    """
    try:
        df = getattr(model, attr_name)
        dict_ = df.loc[obj.id].to_dict()
        obj.notes.update(dict_)
    except KeyError:
        dict_ = {}
        obj.notes = dict_
    except AttributeError:
        pass
    return
