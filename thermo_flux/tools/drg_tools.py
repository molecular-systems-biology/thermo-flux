
#ToDo set up install of equilibrator assests
#from equilibrator_assets.generate_compound import create_compound, get_or_create_compound

from functools import reduce
import pandas as pd
import numpy as np
import numpy.linalg as la

from ..utils.vis import progressbar



from math import floor, log10, inf, copysign


import warnings
def custom_formatwarning(message, category, filename, lineno, line=None):
    return f'{message}\n'
warnings.formatwarning = custom_formatwarning

from cobra.util.array import create_stoichiometric_matrix
from cobra.core import Metabolite
from cobra.core import Reaction as cobra_reaction
from equilibrator_cache.models.compound import Compound as eQ_compound
from equilibrator_cache import PROTON_INCHI_KEY, thermodynamic_constants


import thermo_flux.core.metabolite

from equilibrator_api import Q_, R, FARADAY

from equilibrator_api.component_contribution import find_most_abundant_ms

from equilibrator_api import Reaction


#ToDo set up install of equilibrator assests
def new_eq_metabolite(met):
    """function to add a new metabolite to the local eq database if it does not already exist"""
    met_string = None
    mol_format = None

    if met._unknown:
        return None
    else:
        if hasattr(met.model,'_new_met_index'):
            met.model._new_met_index -= 1
        else:
            met.model._new_met_index = -1

        if 'smiles' in met.annotation:
            met_string = met.annotation['smiles']
            mol_format = 'smiles'

        elif 'InChI' in met.annotation:
            met_string = met.annotation['InChI']
            mol_format = 'inchi'

        if met_string is not None:
            cpd = met.model.lc.get_compounds([met_string], mol_format = mol_format)
            if cpd is not None:
                cpd = cpd[0]
           
        else:
            cpd=None
            met._unknown = True

        return cpd

def get_suitable_ids(met, search = False, update_annotations = False ):
    found = False
    cpd = None
    formula = None
    inchi = None
    searched = False
    annotation = {}

    #if not met.unknown: #still search for metabolites even if unkown... this might slow things down 
       
    # first use metabolite annotations to find metabolite in equilibrator database
    for key, value in met.annotation.items():
        if found:
            break
        if isinstance(value, list):  # if annotations are in a sublist
            for v in value:
                if found:
                    break
                cpd = met.model.cc.get_compound(str(key) + ':' + v)

                if cpd is not None:
                    found = True
                    break
                else: #remove anything after . in identifyer key in case this causes issue 
                    cpd = met.model.cc.get_compound(str(key.split('.')[0]) + ':' + v)
                    if cpd is not None:
                        found = True
                        break

        elif isinstance(value, str):
            cpd = met.model.cc.get_compound(str(key) + ':' + value)
            if cpd is not None:
                found = True
                break
            else: #remove anything after . in identifyer key in case this causes issue 
                    cpd = met.model.cc.get_compound(str(key.split('.')[0]) + ':' + value)
                    if cpd is not None:
                        found = True
                        break
                
    # if the compound has been found but does not have an inchi
    # this suggests it does not have a decomposable structure
    # in this case default to using any InChi key or smiles specified to make a new compoud
    if found:
        if cpd.inchi is None:
            cpd_new = new_eq_metabolite(met)
            if cpd_new is not None:
                found = True
                cpd = cpd_new

    # if the compound was not found at all try and create a new compound based on inchi 
    if not found:
        cpd_new = new_eq_metabolite(met)
        if cpd_new is not None:
            found = True
            cpd = cpd_new

    # as a last resort try searching equilibrator using any common name in metabolite notes. 
    # note this can return incorrect compounds and needs manual checking! 
    if not found:  
        if search:
            if 'common name' in met.notes:
                common_name = met.notes['common name']
                cpd = met.model.cc.search_compound(common_name) 
                searched = True

            if cpd is None and 'charge' not in met.name: ### otherwise a wrong compound can be found for charge metabolite
                cpd = met.model.cc.search_compound(met.name) 
                searched = True

            if cpd is None:
                for key, value in met.annotation.items():
                    if found:
                        break
                    if isinstance(value, list):  # if annotations are in a sublist
                        for v in value:
                            if found:
                                break
                            cpd = met.model.cc.search_compound(value)
                            if cpd is not None:
                                found = True
                                break
                    elif isinstance(value, str):
                        cpd = met.model.cc.search_compound(value)
                        if cpd is not None:
                            found = True
                            break

            if cpd is not None:
                found = True
                        
    if update_annotations: #Warning this can be slow 
        if found:
            annotation = met.annotation
            for identifier in cpd.identifiers:
                namespace = identifier.registry.namespace
                if namespace in annotation:
                    if isinstance(annotation[namespace],list):
                        annotation[namespace].append(identifier.accession)
                    else:
                        annotation[namespace] = [annotation[namespace], identifier.accession]

                else:
                    annotation[identifier.registry.namespace] = identifier.accession

            formula = cpd.formula
            inchi = cpd.inchi

    return cpd, annotation, formula, inchi, searched

def get_compound(met):
    """takes an identifier string, thermo metabolite or equilibrator compound and returns an equilibrator compound"""
    cpd = None
    if type(met) == eQ_compound:
        cpd = met

    if type(met) == str:
        cpd = met.model.cc.get_compound(met)
   
    if type(met) == thermo_flux.core.metabolite.ThermoMetabolite:
        if type(met.compound) == eQ_compound:
            cpd = met.compound
            
        else:
            cpd, annotation, formula, inchi, searched = get_suitable_ids(met)     

    return cpd

def round_and_normalize(numbers, round_dp=2):
    numbers_rounded = np.round(numbers, round_dp)
    sum_rounded = sum(numbers_rounded)
    if sum_rounded < 1:
        max_index = numbers_rounded.argmax()
        numbers_rounded[max_index] += (1 - sum_rounded)
    elif sum_rounded > 1:
        min_index = numbers_rounded.argmin()
        numbers_rounded[min_index] -= (sum_rounded - 1)

    z = [num / sum(numbers_rounded) for num in numbers_rounded]
    return z

def calc_average_charge_protons(compound, pH, pMg, ionic_strength, temperature, accuracy=0.1, round_dp=False,cobra_formula=False):
    """
    Modified from Elad Noor
    Function for calculating microspecies distribution of metabolites
    Microspcies below a certain theshold can be ignored by changing the accuracy. i.e accuracy = 0.1 only considers micorspecies >10% of population 
    """
    cpd = get_compound(compound)
    data = []
    microspecies_df = pd.DataFrame(data, columns=["charge", "number_protons", "number_magnesiums", "ddg_over_rt", "ddg_prime_over_rt"])


    if type(cpd) == eQ_compound:
        if cpd.formula == 'H':
            average_charge = 1
            average_protons = 1
            average_Mg = 0
            return average_charge, average_protons, average_Mg, microspecies_df

        if cpd.inchi_key == PROTON_INCHI_KEY:
            average_charge = 1
            average_protons = 1
            average_Mg = 0
            return average_charge, average_protons, average_Mg, microspecies_df

        if cpd.formula == 'Mg':
            average_charge = 2
            average_protons = 0
            average_Mg = 1
            microspecies_df['number_magnesiums'] = 1
            return average_charge, average_protons, average_Mg, microspecies_df

    if (cpd is None) or (cpd.can_be_transformed() is False) or (cpd.atom_bag is None) or (cpd.atom_bag == {}) or (cobra_formula==True):
        #warnings.warn(f"{compound.id} cannot be transformed due to unknown structure or protons, using COBRA model formula and charge instead", stacklevel=2)
        if type(compound) == thermo_flux.core.metabolite.ThermoMetabolite:
        
            if 'H' in compound.elements:
                average_protons = compound.elements['H']
            else:
                average_protons = 0
            if 'Mg' in compound.elements:
                average_Mg = compound.elements['Mg']
            else:
                average_Mg = 0
            if compound.charge is None:
                average_charge = 0
            else:
                average_charge = compound.charge
            microspecies_df = None
        else:
            average_protons = None
            average_charge = None
            average_Mg = None
            microspecies_df = None

        return average_charge, average_protons, average_Mg, microspecies_df

    
    for ms in cpd.microspecies:
        # perform a Legendre transform (based on Alberty's method) to get the transformed Gibbs energy of each microspecies
        ddg_prime_over_rt = ms.transform(pH=pH.m_as(""), pMg=pMg.m_as(""), ionic_strength_M=ionic_strength.m_as("M"), T_in_K=temperature.m_as("K")
)

        # store all the relevant data in a list
        data.append((ms.charge, ms.number_protons, ms.number_magnesiums, ms.ddg_over_rt, ddg_prime_over_rt))

    # convert the list of values into a DataFrame
    microspecies_df = pd.DataFrame(data, columns=["charge", "number_protons", "number_magnesiums", "ddg_over_rt", "ddg_prime_over_rt"])

    # mark the row which corresponds to the major MS at pH 7 (according to ChemAxon)
    # this microspecies should always have a ddg_over_rt equal to 0, but after the transform it might not be 0 anymore
    #microspecies_df["is_major_ms"] = (microspecies_df.charge == cpd.net_charge) & (microspecies_df.number_magnesiums == 0) #depreciated at equilibrator v6?

    if len(microspecies_df) == 0:
        average_charge = None
        average_protons = None
        average_Mg = None
        return average_charge, average_protons, average_Mg, microspecies_df
    
    if cpd.formula == 'H':
        average_charge = 1
        average_protons = 1
        average_Mg = 0
        return average_charge, average_protons, average_Mg, microspecies_df

    if cpd.inchi_key == PROTON_INCHI_KEY:
        average_charge = 1
        average_protons = 1
        average_Mg = 0
        return average_charge, average_protons, average_Mg, microspecies_df

    if cpd.formula == 'Mg':
        average_charge = 2
        average_protons = 0
        average_Mg = 1
        microspecies_df['number_magnesiums'] = 1
        return average_charge, average_protons, average_Mg, microspecies_df
    
    # the abundances are proportional to e^(-ΔG'), so all we need is to calculate all the exponents and then normalize the sum to 1.
    # however, this is numerically unstable because the log-values can span a very wide range. therefore we normalize the values
    # in log-space using the logaddexp function, and only then calculate the exponents

    # calculate the normalizing factor for the Boltzmann distribution (log of the sum of exponents)
    log_boltzmann_normalizing_factor = reduce(np.logaddexp, -microspecies_df.ddg_prime_over_rt.values)

    microspecies_df["abundance"] = np.exp(-microspecies_df.ddg_prime_over_rt.values - log_boltzmann_normalizing_factor)

    #drop microspecies < 10% of distribution and renormalise to 1.0
    abundance_round = []
    for abundance in (microspecies_df["abundance"]):
        if abundance < accuracy:
            abundance_round.append(0) 
        else: 
            abundance_round.append(abundance)
            
    #renormalise to 1.0 
    abundance_round_norm = [(float(i)/sum(abundance_round)) for i in abundance_round] 
         
    microspecies_df["abundance_norm"] = abundance_round_norm

    if round_dp is not False:
        microspecies_df["abundance_round"] = round_and_normalize(abundance_round_norm, round_dp=round_dp) 
    else:
        microspecies_df["abundance_round"] = abundance_round_norm

    average_charge = sum(microspecies_df["abundance_round"]*microspecies_df["charge"])
    average_protons = sum(microspecies_df["abundance_round"]*microspecies_df["number_protons"])
    average_Mg = sum(microspecies_df["abundance_round"]*microspecies_df["number_magnesiums"])

    if round_dp is not False:
        average_charge = round(average_charge, round_dp)
        average_protons = round(average_protons, round_dp)
        average_Mg = round(average_Mg, round_dp)

    return average_charge, average_protons, average_Mg, microspecies_df

def major_microspecies(met, pH, pMg, ionic_strength_M, T_in_K ,cobra_formula=False):
    """Calculates the major microspecies of a compound. Useful for transporter calculations"""

    cpd = get_compound(met)

    if cpd is None or  cpd.can_be_transformed() == False or cobra_formula==True:
        #ToDo make this a single function to return COBRA model information 
       # warnings.warn(f"{met.id} cannot be transformed due to unknown structure or protons, using COBRA model formula and charge instead", stacklevel=2)
        if type(met) == thermo_flux.core.metabolite.ThermoMetabolite:
            if 'H' in met.elements:
                protons = met.elements['H']
            else:
                protons = 0
            if 'Mg' in met.elements:
                Mg = met.elements['Mg']
            else:
                Mg = 0
            if met.charge is None:
                charge = 0
            else:
                charge = met.charge

        else:
            protons = None
            charge = None
            Mg = None
        return charge, protons, Mg 


    elif cpd.inchi_key == PROTON_INCHI_KEY:
        protons = 1
        charge = 1
        Mg=0
    else:
        ms = find_most_abundant_ms(cpd, pH, pMg, ionic_strength_M, T_in_K)
        protons = ms.number_protons
        charge = ms.charge
        Mg = ms.number_magnesiums 
    return charge, protons, Mg 

def pka_graph(metabolite, pMg = None, ionic_strength = None, temperature = None, accuracy=0, round_dp=False):
    '''returns a dataframe of the charge distribution of a metabolite at different pHs
    Parameters
    ----------
    metabolite: cobra metabolite object
    pMg: float or pint quantity
    ionic_strength: pint quantity

    Returns
    -------
    dataframe with pH as index and charge as columns
    '''
  
    if pMg is None:
        pMg = metabolite.model.pMg[metabolite.compartment]

    if ionic_strength is None:
        ionic_strength = metabolite.model.I[metabolite.compartment]

    if temperature is None:
        temperature = metabolite.model.T


    ph_data = {}
    for pH in np.linspace(0,14,100):
        pH = Q_(pH,'')
        charge, protons, mg, df = thermo_flux.tools.drg_tools.calc_average_charge_protons(metabolite, pH=pH, pMg = pMg, ionic_strength = ionic_strength, temperature = temperature, accuracy=accuracy, round_dp=round_dp)
        df['charge_protons_magnesiums'] = list(zip(df['charge'], df['number_protons'], df['number_magnesiums']))
        df['charge_protons_magnesiums'] = list(zip(df['charge'], df['number_protons'], df['number_magnesiums']))
        df = df[['abundance_round','charge_protons_magnesiums']].set_index('charge_protons_magnesiums')

        ph_data[pH.m] = df

    df = pd.concat(ph_data, axis=1)

    df.columns = df.columns.droplevel(1)
    df = df.T

    #drop columns where all values <0.01
    df = df.loc[:, (df > 0.01).any(axis=0)]

    return df

def calc_dfG0(tmodel, fit_unknown_dfG0=False):
    default_dfg0_std = tmodel.rmse_inf

    known_met_idx = []
    unknown_met_idx = []
    unknown_mets = []

    dfG0_mean = Q_(np.zeros((len(tmodel.metabolites), 1)),"kJ/mol")
    dfG0_cov_sqrt = Q_(np.zeros((len(tmodel.metabolites), len(tmodel.metabolites))),"kJ/mol")


    known_met_list = []

    for i, met in enumerate(tmodel.metabolites):
        if met.compound is None:
            unknown_met_idx.append(i)
            unknown_mets.append(met)
        else:
            known_met_idx.append(i)
            known_met_list.append(Reaction({met.compound:1}))

    #calculate dfG of metabolites known to equilibrator 
    dfG0_mean_known, dfG0_cov_sqrt_known = tmodel.cc.standard_dg_multi(known_met_list,"fullrank")

    dfG0_mean[known_met_idx,0] = dfG0_mean_known
    dfG0_cov_sqrt[known_met_idx, 0:dfG0_cov_sqrt_known.shape[1]] = dfG0_cov_sqrt_known


    #set unknown metabolites to have a reasonable formation energy
    #is the mean of other formation energies reasonable? 
    #would 0 be a better estimate?
    #is this very case specific i.e generic proteins etc. are not normal metabolites 
    default_dfg0_mean = np.mean(dfG0_mean[known_met_idx]) #this value is not needed becuase of the fitting later?
    dfG0_mean[unknown_met_idx] = Q_(0,"kJ/mol") #set unknown metabolites to have formation energy of 0 unless specified in tmet
        
    
    #for metabolites which have an identifier in the equilibrator cache but we want to ignore it and change the uncertainty
    #only for redox carrier and biomass can we set our own formation energy 
    #aslo for photons 
    #can only do this for metabolites that are not correlated with any others (apart from themselves)
    
    special_mets_idx = []    
    for i, met in enumerate(tmodel.metabolites):
        if any([met.redox, met.biomass,]):
            if len(np.nonzero(dfG0_cov_sqrt[i])[0]) <= 1:   #only metabolites with a single entry in dfG0_cov_sqrt
                dfG_mean = Q_(0,"kJ/mol")
                if met.dfG0 is not None: 
                    dfG_mean = met.dfG0 #this shold be a feature of a reaction not a metabolite - stoichometry need accounting for?
                
                #if only a value for dfG0prime has been set then use this for the dfG0 value otherwise dfG_mean should represent untransformed value
                if all([(met.dfG0prime is not None), (met.dfG0 is None)]):  
                    dfG_mean = met.dfG0prime 

                special_mets_idx.append(i)
                dfG0_mean[i] = dfG_mean
                
                if met.dfG_SE is None:
                    dfG_SE = Q_(0,"kJ/mol")
                    print(met.id, dfG_SE)
                else:
                    dfG_SE = met.dfG_SE

                unknown_basis_idx = (np.abs(dfG0_cov_sqrt[i])).argmax()
                dfG0_cov_sqrt[i,unknown_basis_idx] = dfG_SE
                
    #identify metabolites for whcih equilibrator has no estimate of dfG 
    unknown_basis_idx = [i for i in range(dfG0_cov_sqrt.shape[1]) if np.any(np.abs(dfG0_cov_sqrt[:, i]) >= default_dfg0_std)]
    
    #enforce correlation of multiple instances of unknown metabolites in different compartments 
    dfG_basis_size = dfG0_cov_sqrt_known.shape[1]
    unknown_met_basis_id = {}
    

    for met_idx in unknown_met_idx:
        if met_idx not in special_mets_idx: #ignore special metabolites that have defined dfGs 
            met = tmodel.metabolites[met_idx]
            if met.accession in unknown_met_basis_id: #if the metabolite already has column in the unknown basis
                if met.dfG_SE is None:
                    dfG_SE = default_dfg0_std
                else:
                    dfG_SE = default_dfg0_std 

                dfG0_cov_sqrt[met_idx, unknown_met_basis_id[met.accession]] = dfG_SE

            else:#otherwise add a new basis for the unknown metabolite 

                unknown_met_basis_id[met.accession] = dfG_basis_size
                if met.dfG_SE is None:
                    dfG_SE = default_dfg0_std
                else:
                    dfG_SE = met.dfG_SE
                dfG0_cov_sqrt[met_idx, dfG_basis_size] = dfG_SE
                dfG_basis_size += 1


    #return unknown basis for fitting later
    unknown_basis_idx += range(dfG0_cov_sqrt_known.shape[1], dfG_basis_size)
    unknowns_basis = dfG0_cov_sqrt[:, unknown_basis_idx]
    unknown_shift = Q_(0, 'kJ/mol')

    if fit_unknown_dfG0:
        S = create_stoichiometric_matrix(tmodel)

        unknown_shift = _fit_unknown_dfG(S, unknowns_basis, dfG0_mean[:, 0])

    dfG0_mean = dfG0_mean[:, 0]+unknown_shift

    #drop degrees of freedom that are all 0 
    dfG0_cov_sqrt = dfG0_cov_sqrt[:,np.where(dfG0_cov_sqrt.any(axis=0))[0]]

    return dfG0_mean, dfG0_cov_sqrt, unknowns_basis

def _fit_unknown_dfG(S, unknowns_basis, dfG0_prime):
    X = S.T @ unknowns_basis.m  # this is unknown part of formation energies 
    y = -S.T @ dfG0_prime.m_as("kJ/mol")  # these are drG0prime

    # try to minimise
    beta = la.lstsq(X, y, rcond=None)[0]  #shift dfG of completely unknown metabolites to minmise effect of unknown dfG on drG estimation?
    return Q_(unknowns_basis @ beta, "kJ/mol")


def calc_dfG_transform(met):
    #calucalte transform required for every metabolite based on model compartment information 
    cpd = get_compound(met)

    pH = met.model.pH[met.compartment]
    I = met.model.I[met.compartment]
    T = met.model.T
    pMg = met.model.pMg[met.compartment]
    
    if cpd is None:
        if met.elements != {} and 'H' in met.elements.keys(): #if the metabolite has a cobra formula/charge then it can be transformed witht the number of protons
            charge = 0 #default to 0 charge unless otherwise specified in the met.charge 
            if met.charge:
                charge = met.charge
            dfG_transform = thermodynamic_constants._legendre_transform(
            pH= pH.m,
            pMg= 14,#ignore magnesium 
            ionic_strength_M= I.m_as('M'),
            T_in_K= T.m,
            charge= charge,
            num_protons= met.elements['H'], 
            num_magnesiums= 0
            )*(R*T)
        else : 
            dfG_transform = Q_(0,"kJ/mol")
    elif cpd.can_be_transformed():
        dfG_transform = cpd.transform(pH, I, T, pMg)
    else:
        if met.elements != {} and 'H' in met.elements.keys(): #if the metabolite has a cobra formula/charge then it can be transformed with the number of protons
            charge = 0 #default to 0 charge unless otherwise specified in the met.charge 
            if met.charge:
                charge = met.charge
            dfG_transform = thermodynamic_constants._legendre_transform(
            pH= pH.m,
            pMg= 14,#ignore magnesium 
            ionic_strength_M= I.m_as('M'),
            T_in_K= T.m,
            charge= charge,
            num_protons= met.elements['H'], 
            num_magnesiums= 0
            )*(R*T)
        else : 
            dfG_transform = Q_(0,"kJ/mol")

      #for biomass manually transform based on the formula 
    if met.biomass:
        dfG_transform = thermodynamic_constants._legendre_transform(
            pH= pH.m,
            pMg= 14,#ignore magnesium 
            ionic_strength_M= I.m_as('M'),
            T_in_K= T.m,
            charge= 0,
            num_protons= met.average_charge_protons()[1], 
            num_magnesiums= 0
        )*(R*T)

    return dfG_transform
       

def calc_dfG0prime(tmodel, fit_unknown_dfG0=False):
    dfG0_mean, dfG0_cov_sqrt, unknowns_basis = calc_dfG0(tmodel, fit_unknown_dfG0=fit_unknown_dfG0)
    
    dfG_transforms = []
    for met in progressbar(tmodel.metabolites, '', 40, item_label_attribute = 'id'):
        dfG_transform = calc_dfG_transform(met)
        dfG_transforms.append(dfG_transform.m)

    dfG0prime_mean = dfG0_mean + Q_(np.array(dfG_transforms), "kJ/mol")

    return dfG0_mean, dfG0prime_mean, dfG0_cov_sqrt, unknowns_basis
    

def calc_drG0(S, dfG0):
    drG0 = S.T @ dfG0
    return drG0

def calc_model_drG0(tmodel):
    S = create_stoichiometric_matrix(tmodel)
    dfG0 = Q_(np.array([met.dfG0.m for met in tmodel.metabolites]), 'kJ/mol')
    drG0 = calc_drG0(S, dfG0)
    return drG0

def calc_model_drG0prime(tmodel):
    dfG0prime = Q_(np.array([met.dfG0prime.m for met in tmodel.metabolites]), 'kJ/mol')
    S = create_stoichiometric_matrix(tmodel)
    drG0prime = calc_drG0(S, dfG0prime)

    return drG0prime

def calc_phys_correction(tmodel):
    S = create_stoichiometric_matrix(tmodel)
    phys_conc = np.array([np.log(0.001) if not met.ignore_conc else 0 for met in tmodel.metabolites])
    phys_conc_correction = R * tmodel.T * S.T @ phys_conc

    return phys_conc_correction

def formula_dict_to_string(formula):
    formula_string = ''
    for key, value in formula.items():
        if key not in ['charge', 'e-']:
            formula_string += key + str(value)
    return formula_string

def calc_biomass_formula(biomass_rxn):
    '''calculate biomass elemental compositon from biomass equation'''

    #reset biomass formula to be empty
    for met in biomass_rxn.metabolites:
        if met.biomass:
            met.formula = ''
            met.elements = {}
            comp = met.compartment

    biomass_atom_bag = {}
    for met, stoich in biomass_rxn.metabolites.items():
        if met.compound is not None:
            for element, value in met.compound.atom_bag.items():
                if element in biomass_atom_bag:
                    biomass_atom_bag[element] += -1*stoich*value
                else:
                    biomass_atom_bag[element] = -1*stoich*value
        else:
            for element, value in met.elements.items():
                if element in biomass_atom_bag:
                    biomass_atom_bag[element] += -1*stoich*value
                else:
                    biomass_atom_bag[element] = -1*stoich*value

        #should net_elements be used here to calcualte protons and charge in biomass? accounting for the species of the precursors?
        biomass_protons = -1*biomass_rxn.net_elements()[0][biomass_rxn.model.proton_dict[comp]]
        biomass_atom_bag['H'] = biomass_protons
            
        #round formula to help with long floating point numbers
        for element in biomass_atom_bag:
            biomass_atom_bag[element] = np.round(biomass_atom_bag[element], 10)

    biomass_formula_string = formula_dict_to_string(biomass_atom_bag)

    return biomass_formula_string, biomass_atom_bag

def calculate_biomass_dfG0(biomass):
    '''calculate the formation energy of biomass from the biomass formula
    
    Parameters
    ----------
    biomass: cobra metabolite object

    Returns
    -------
    dfG0_bm: Quantity
        formation energy of biomass units are defined as kJ/mol for compatibility with other reactions
        however, actual units for this value are J/gDW due to conversion in biomass eqation from mmol to gDW
    '''

    if 'C' in biomass.elements.keys(): 
        if biomass.elements['C'] > 0:
            biomass_formula = biomass.elements
    else:
        #idenitfy biomass reaction from biomass metabolite
        #assume the reaction biomass metabolites is involved in with the most metabolites is the biomass reaction 
        biomass_rxn = [rxn for rxn in biomass.reactions][[len(rxn.metabolites) for rxn  in biomass.reactions].index(max([len(rxn.metabolites) for rxn  in biomass.reactions]))]
        biomass_formula_string, biomass_formula = calc_biomass_formula(biomass_rxn)
        biomass.formula = biomass_formula_string
        
    #normalise to cmol 
    cmol = biomass_formula['C']
    for element in biomass_formula:
        biomass_formula[element] = biomass_formula[element]/cmol

    #assume biomass is not charged
    biomass_formula['charge'] = 0

    #estimate biomass formation energy from biomass formula using Battley 1993 method
    dfGbm_g = thermo_flux.tools.drg_tools.dfGbm(biomass_formula, units = 'kJ/g')[0]

    #for unit parity with metabolic reactions we have convert to kJ/mol although these are not the actual units
    #fluxes are in mmol/gDW/h
    #biomass reaction is in units of gDW/gDW/h
    #biomass Gdiss is in units of J/gDW/h

    #the transform is already in kJ/mol and depends on the number of H in one mol of biomass
    #since 1 mol H weighs 1 g, the transform is equivalent to kJ/gDW

    dfG0_bm = Q_(dfGbm_g.m*1000, 'kJ/mol') # *1000 to convert to J as fluxes are in mmol/gDW/h

    return dfG0_bm

def dfGbm(formula = {}, units = 'kJ/mol', Mr_bio = None):
    '''
      Calculates the formation energy of biomass or macromolecules based on their empirical formula
      modified from Saadat et. al Entropy 2020, 22(3), 277; https://doi.org/10.3390/e22030277 https://gitlab.com/qtb-hhu/thermodynamics-in-genome-scale-models/-/blob/master/EnergyOfFormationBiomass.py?ref_type=heads

      Parameters
      ----------
      formula : dict
            Empirical formula of macromolecule.
      units : str, optional
            Units of the output. The default is 'kJ/mol'.
      Mr_bio : Quantity, optional
            Molecular weight of biomass in units carbon mol/gDW . The default 
            is None and will be automatically calcualted from the 
            elemental composition.
      Returns
      -------
        Gf : Quantity
            Gibbs energy of formation in units of 'units'      
        Gc : Quantity
            Gibbs energy of combustion in units of 'units'
        y : float
            Degree of reduction
        stoich : array
            Stoichiometry of combustion reaction
    '''

    C = formula['C'] if 'C' in formula else 1
    H = formula['H'] if 'H' in formula else 0
    O = formula['O'] if 'O' in formula else 0
    N = formula['N'] if 'N' in formula else 0
    P = formula['P'] if 'P' in formula else 0
    S = formula['S'] if 'S' in formula else 0
    K = formula['K'] if 'K' in formula else 0
    Mg = formula['Mg'] if 'Mg' in formula else 0
    Ca = formula['Ca'] if 'Ca' in formula else 0
    Fe = formula['Fe'] if 'Fe' in formula else 0
    charge = formula['charge'] if 'charge' in formula else 0

    #Energies of formation for various chemical compounds (see Battley 1993)
    Gf_B = {'CO2_g':-394.36,
            'N2_g': 0,
            'O2_g':0,
            'KOH_c':-379.11,
            'K2SO4_c':-1321.43,
            'P4O10_c':-2697.84,
            'H2O_lq':-237.18}

    # Molecular weight of biomass in g/mol
    if 'mol' in units:
        Mr_bio = 1
    elif units == 'kJ/g':
        if Mr_bio is not None: #if user defined denisty (Cmol/g)
            pass
        else:
            Mr_bio = (C*12.011 + H*1.008 + O*15.999 + N*14.007 + P*31.0 + S*32.06 + K*39.098 + Mg*24.305 + Ca*40.078 + Fe*55.845)

    #Assumes combustion reaction see Battley 1993
    #biomass + a O2 + b KOH -> c CO2 + d N2 + e P4O10 + f K2SO4 + g H20   
    a = np.array([
        [0,0,1,0,0,0,0], #C
        [0,1,0,0,0,0,2], #H
        [2,1,2,0,10,4,1], #O
        [0,0,0,2,0,0,0], #N
        [0,0,0,0,4,0,0], #P
        [0,0,0,0,0,1,0], #S
        [0,1,0,0,0,2,0]]) #K

    b = [C,H,O,N,P,S,K]
    stoich = np.linalg.solve(a,b)

    #calculate degree of reduction 
    y = 4*C + H - 2*O + 0*N + 5*P + 6*S - charge 
        
    # Calculates energy of combustion according to Battley1993. 
    # Other formulas can be found in Stockar1993 and Minkevich1973 
    Gc = y*-107.90

    Gf = (Gc - stoich[0]*Gf_B['O2_g'] - stoich[1]*Gf_B['KOH_c'] - 
                stoich[2]*Gf_B['CO2_g'] - stoich[3]*Gf_B['N2_g'] - 
                stoich[4]*Gf_B['P4O10_c'] - stoich[5]*Gf_B['K2SO4_c'] - 
                stoich[6]*Gf_B['H2O_lq'])

    return Q_(-Gf/Mr_bio, units),Q_(Gc/Mr_bio, units),y,stoich #minus Gf because biomass is a substrate in the combustion reaction         

def proton_dict(tmodel):
    """identifies protons in the model so they can be added to reactions for balancing """

    proton_annotations = {'Kegg': 'C00080',
     'bigg.metabolite': 'h',
     'chebi': 'CHEBI:15378',
     'hmdb': 'HMDB59597',
     'kegg.compound': 'C00080',
     'metacyc.compound': 'PROTON',
     'reactome': '1132304',
     'sabiork.compound': '39',
     'seed': 'cpd00067',
     'metanetx.chemical': 'MNXM1',
     'synonyms': 'H+'}

    proton_dict = {}
    compartments = []
    for met in tmodel.metabolites:
        for key, value in proton_annotations.items():
            if key in met.annotation:
                if met.annotation[key] == value:
                    proton_dict[met.compartment] = met
                    compartments.append(met.compartment)
                    pass
                #if annotations are in list form 
                if met.annotation[key] is list:
                    if value in met.annotation[key]:
                        proton_dict[met.compartment] = met
                        compartments.append(met.compartment)
                        pass

    #ToDo add case for when there are no protons in the model at all
    if len(proton_dict) == 0:
        for comp in tmodel.compartments: 
            proton = thermo_flux.core.metabolite.ThermoMetabolite(Metabolite('h_'+comp, 'H', 'h',1,comp))
            proton.annotation = proton_annotations
            tmodel.add_metabolites([proton])
            proton_dict[comp] = proton
    else:
        #if some protons are in the model but not all compartments add new metabolites 
        for comp in tmodel.compartments:
            if comp not in proton_dict:
                proton = (list(proton_dict.values())[0]).copy() 
                proton.compartment = comp
                proton.id = proton.id[:-2]+'_'+comp
                proton.compound = None
                proton._model = tmodel
                proton.compound
                tmodel.add_metabolites([proton])
                proton_dict[comp] = proton

    return proton_dict

def charge_dict(tmodel):
    """identifies charge in the model so they can be added to reactions for balancing """
    charge_dict = {}
    compartments = []
    for met in tmodel.metabolites:
        if 'charge' in met.id:
            charge_dict[met.compartment] = met
            compartments.append(met.compartment)
            met.ignore_conc = True

        
    #ToDo add case for when there is no charge in the model at all
    if len(charge_dict) == 0:
        for comp in tmodel.compartments: 
            charge = thermo_flux.core.metabolite.ThermoMetabolite(Metabolite('charge_'+comp, '', 'charge',1,comp))
            tmodel.add_metabolites([charge])
            charge_dict[comp] = charge
            charge.ignore_conc = True


    else: 
        #if some charges are in the model but not all compartments add new metabolites 
        for comp in tmodel.compartments:
            if comp not in charge_dict:
                charge = (list(charge_dict.values())[0]).copy() 
                charge.compartment = comp
                charge.id = charge.id[:-2]+'_'+comp
                charge_dict[comp] = charge
                charge.ignore_conc = True

    return charge_dict

def mg_dict(tmodel):
    """identifies magnesiums in the model so they can be added to reactions for balancing """

    mg_annotations = {'bigg.metabolite': 'mg2',
         'chebi': 'CHEBI:18420',
         'hmdb': 'HMDB00547',
         'kegg': 'C00305',
         'metacyc.compound': 'MG+2',
         'reactome': '109496',
         'sabiork.compound': '1327',
         'seed': 'cpd00254',
         'metanetx.chemical': 'MNXM106383',
         'synonyms': 'Magnesium'}

    mg_dict = {}
    compartments = []
    for met in tmodel.metabolites:
        if met.formula == 'Mg1' or met.formula == 'Mg':
            mg_dict[met.compartment] = met
            compartments.append(met.compartment)
        else:
            for key, value in mg_annotations.items():
                if key in met.annotation:
                    if met.annotation[key] == value:
                        mg_dict[met.compartment] = met
                        compartments.append(met.compartment)
                        pass
                     #if annotations are in list form 
                    if met.annotation[key] is list:
                        if value in met.annotation[key]:
                            mg_dict[met.compartment] = met
                            compartments.append(met.compartment)
                            pass

    if len(mg_dict) == 0:
        for comp in tmodel.compartments: 
            mg = thermo_flux.core.metabolite.ThermoMetabolite(Metabolite('Mg_'+comp, 'Mg', 'Mg',2,comp))
            mg.annotation = mg_annotations
            tmodel.add_metabolites([mg])
            mg_dict[comp] = mg
 
    else:
        #if some mg are in the model but not all compartments add new metabolites 
        for comp in tmodel.compartments:
            if comp not in mg_dict:
                mg = (list(mg_dict.values())[0]).copy() 
                mg.compartment = comp
                mg.id = mg.id[:-2]+'_'+comp
                mg_dict[comp] = mg

    return mg_dict

def _transport_direction(reaction):
    """Function to define the inner and outer compartments of a transport reaction. 
    Useful for consistently defining the major microspecies that is transported considering the conditions 
    of the inner compartment"""

    phi_dict = reaction.model.phi_dict

    #reaction = reaction.copy()

    compartments = list(reaction.compartments)

    if len(compartments)>1:
        inner_comp = reaction.reactants[0].compartment 
        compartments.remove(inner_comp)
        outer_comp = compartments[0]
    else:
        inner_comp = compartments[0]
        outer_comp = compartments[0]

    # Convention for membrane potentials being V_inside - V_outside
    # potential is generally higher outside so if potential is -'ve then inner and outer are correct 
    # if potential is +'ve then inner and outer should be swapped 
    e_potential_difference = phi_dict[inner_comp][outer_comp]
    swapped = 1
    if e_potential_difference < 0:
        inner_comp_temp = inner_comp
        inner_comp = outer_comp
        outer_comp = inner_comp_temp

        e_potential_difference = phi_dict[inner_comp][outer_comp]
        metabolites = reaction.metabolites
        swapped = {}
        for met, stoich in metabolites.items():
            swapped[met] = stoich*-1
            
        #reaction.add_metabolites(swapped)
        #reaction.add_metabolites(swapped)
        swapped = -1

    if len(reaction.compartments) == 1:
        outer_comp = None

    return reaction, swapped, inner_comp, outer_comp

#separate reaction into inner and outer half reactions 
def comp_split(reaction, compartment):
    metabolites = {}
    for metabolite in reaction.metabolites:
        if metabolite.compartment == compartment:
            metabolites[metabolite] = reaction.metabolites[metabolite]
    return (metabolites)


def transported_c_h(reaction, round_dp = False, verbose=False, rxn_already_balanced = True):
    reaction_swapped, swapped, inner_comp, outer_comp = _transport_direction(reaction)
    transported_mets = calc_transported_mets(reaction)
    

    #define compartment conditions
    pMg = Q_(14,'') #for transported mets set low Mg conc as we assume Mg is not transported 
    pH_inner = reaction.model.pH[inner_comp]
    ionic_strength_inner = reaction.model.I[inner_comp]
    #pH_outer = reaction.model.pH[outer_comp]
    #ionic_strength_outer = reaction.model.I[outer_comp]
    temperature = reaction.model.T

    inner_mets = comp_split(reaction,inner_comp)
    outer_mets = comp_split(reaction,outer_comp)

    #calculate net protons and charge of inner and outer metabolites and check for balance 
    #for protons and charge multiply by stoichiometry 
    n_h_inner = 0
    z_inner = 0
    e_inner = 0
    n_h_outer = 0
    z_outer = 0
    e_outer = 0 

    # specific case for transporters that involve chemical transformation 
    # user defined transported metabolites overrides automatic calculation
    if reaction.transported_mets is not None:
        
        for met, stoich in reaction.transported_mets.items():
            charge, protons, Mg = major_microspecies(met, pH_inner, pMg, ionic_strength_inner, temperature,cobra_formula=reaction._cobra_formula)

            if met.compartment == inner_comp:
                swapped = 1
            elif met.compartment == outer_comp:
                swapped = -1

            n_h_inner += swapped*stoich * protons
            z_inner +=  swapped*stoich * charge

            n_h_outer -= swapped*stoich * protons
            z_outer -=  swapped*stoich * charge

        balanced = True
        return n_h_inner, z_inner,n_h_outer, z_outer,inner_comp, outer_comp, balanced 
   
    #first check if any metabolites don't have a equilibrator formula and need to get cobra formula
    ##then we have to use cobra formulas for all compounds in the reaction
    for met, stoich in inner_mets.items(): 
        if met not in reaction.model.charge_dict.values():
             
            #if met in reaction.model.proton_dict().values():
            #   if reaction.transported_h is not None:
            #       stoich = reaction.transported_h[inner_comp]
            if not rxn_already_balanced:
                if met in transported_mets:
                    charge, protons, Mg = major_microspecies(met, pH_inner, pMg, ionic_strength_inner, temperature)
                else:
                    charge, protons, Mg, ms_df = met.average_charge_protons(pMg=pMg, round_dp = round_dp) #major ms in metabolite compartment ignore mg ToDo check this doesn't cause gibbs enery balance issues if balancing magneisum... 
            else:
                #ToDo this ovverides major microespeceis calcualtion - makes this func work with balanced reactions
                charge, protons, Mg, ms_df = met.average_charge_protons(pMg=pMg, round_dp = round_dp) #major ms in metabolite compartment ignore mg ToDo check this doesn't cause gibbs enery balance issues if balancing magneisum... 

            #calculate net electrons and add to charge term 
         #   e = 0 
          #  if met.compartment is not None:
           #     if met.compound.atom_bag is not None:
            #        if 'e-' in met.compound.atom_bag:
             #           e = met.compound.atom_bag['e-']  

            #z_inner -= stoich*e
        #    print(met.id, stoich *protons,n_h_inner)
            n_h_inner += stoich * protons
            z_inner +=  stoich * charge
          #  print(met.id, stoich *protons,n_h_inner)
        # if cobra_formula:
        #     print('using cobra formula',met,charge,protons)
    for met, stoich in outer_mets.items():
        if met not in reaction.model.charge_dict.values():

            #if met in reaction.model.proton_dict().values():
            #   if reaction.transported_h is not None:
            #       stoich = reaction.transported_h[outer_comp]

            if not rxn_already_balanced:

                if met in transported_mets:
                    charge, protons, Mg = major_microspecies(met, pH_inner, pMg, ionic_strength_inner, temperature) #consider inner conditions for transported metabolites 
                else:
                    charge, protons, Mg, ms_df = met.average_charge_protons(pMg=pMg, round_dp = round_dp) #major ms in metabolite compartment ignore mg
            else: 
            #ToDo this ovverides major microespeceis calcualtion - makes this func work with balanced reactions
                charge, protons, Mg, ms_df = met.average_charge_protons(pMg=pMg, round_dp = round_dp) #major ms in metabolite compartment ignore mg


            #calculate net electrons and add to charge term 
           # e = 0 
            #if met.compartment is not None:
             #   if met.compound.atom_bag is not None:
              #      if 'e-' in met.compound.atom_bag:
               #         e = met.compound.atom_bag['e-']
                  
            #z_outer -= stoich*e

          #  print(met.id, stoich *protons,n_h_outer)

            n_h_outer += stoich * protons
            z_outer +=  stoich * charge
          #  print(met.id, stoich *protons,n_h_outer)

    # special case for transporters involving ion transport 
    # free ions are defined in the reaction._transported_charge dict 
    if reaction.transported_charge is not None:
            free_charge_trans_inner = reaction._transported_charge[inner_comp]
            free_charge_trans_outer = reaction._transported_charge[outer_comp]

            z_inner += free_charge_trans_inner#add new transported free charge  
            z_outer += free_charge_trans_outer 


    tolerance = 1e-6

    if (abs(n_h_inner + n_h_outer) > tolerance) or (abs(z_inner + z_outer) > tolerance):
        balanced = False
        balance = "unbalanced"
    else:
        balanced = True
        balance = "balanced"

    #if balanced == False:
    #   warnings.warn("inner and outer net half-reactions of "
    #                   f"{reaction.id} are "
    #                   f"{balance}:"
    #                   f"n_h(inner) = {n_h_inner}, n_h(outer) = {n_h_outer}, "
    #                   f"z(inner) = {z_inner}, z(outer) = {z_outer}, ", stacklevel=2)

    reaction.balanced=balanced ## add balanced flag to retrieve later for reporting

    return n_h_inner, z_inner, n_h_outer, z_outer,inner_comp, outer_comp, balanced

def calc_transported_mets(reaction):
    """calculates metabolites transported in a transport reaction"""
    if len(reaction.compartments) != 2:
        return {}
    
    reaction_swapped, swapped, inner_comp, outer_comp = _transport_direction(reaction)

    inner_mets = comp_split(reaction,inner_comp)
    outer_mets = comp_split(reaction,outer_comp)
   
    transported_mets = {}
    
    #define metabolites transported as those that appear in both compartments
    for met, stoich in inner_mets.items():
        if met.accession in [met_out.accession for met_out in outer_mets]: 
            stoich_in = stoich
            for met_out, stoich in outer_mets.items():
                if met.accession == met_out.accession:
                    stoich_out = stoich
                    #if charge or protons or mg then take stoichiometry 
                    if any([met in reaction.model.proton_dict.values(),
                            met in reaction.model.charge_dict.values(),
                            met in reaction.model.mg_dict.values()]):
                        transported_mets[met] = stoich_in
                        transported_mets[met_out] = stoich_out
                    else: # take min value and keep sign 
                        transported_mets[met] = np.sign(stoich_in)*min([abs(stoich_in), abs(stoich_out)])
                        transported_mets[met_out] = np.sign(stoich_out)*min([abs(stoich_in), abs(stoich_out)])


    if reaction.transported_mets is not None:
        for met, stoich in reaction.transported_mets.items():
            transported_mets[met] = stoich

                            
    return transported_mets


def calc_drGtransport(reaction, round_dp = False, rxn_already_balanced=True):
    if len(reaction.compartments) == 2:
        transported_protons, transported_charge, outer_h,outer_z,inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_dp,rxn_already_balanced=rxn_already_balanced)
        #proton stoichiomtry
        proton_stoich = {}
        for met, stoich in reaction.metabolites.items():
            if met in reaction.model.proton_dict.values():
                proton_stoich[met] = stoich
        #charge stoichiometry
        charge_stoich = {}
        for met, stoich in reaction.metabolites.items():
            if met in reaction.model.charge_dict.values():
                charge_stoich[met] = stoich

        #transport reactions must be balanced to accurately caclulate drGtransport 
        if (balanced is False): 
            reaction_balance(reaction, balance_charge = False, balance_mg = False, round_dp = round_dp,rxn_already_balanced=False)
            transported_protons, transported_charge,outer_h,outer_z,inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_dp,rxn_already_balanced=True)
            warnings.warn(f"WARNING: {reaction.id} was automatically balanced to allow accurate calculation of drG transport. The new reaction stoichiometry is: {reaction.reaction}. \n Total transported H: {transported_protons}, total transported charge: {transported_charge}, Transported free H+: {reaction.transported_h}", stacklevel=2)

        #if balanced is still false after this raise a warning that the reaction could not be balanced and needs manually curating
        #this occurs if the imported model is incorectly balanced or the pH of the imported model does not match the pH of the thermo model being updated
        if not balanced:
            warnings.warn(f"WARNING: {reaction.id} is not balanced and could not be automatically balanced, please check reaction stoichiometry", stacklevel=2)
            #return stoichiometry back to original imported model
            #remove any protons and charge from reaction and replace with original stoichiometry
            for met, stoich in reaction.metabolites.items():
                if met in reaction.model.proton_dict.values():
                    reaction.add_metabolites({met:-stoich})
                if met  in reaction.model.charge_dict.values():
                    reaction.add_metabolites({met:-stoich})
            #add back in original proton and charge stoichiometry
            reaction.add_metabolites(proton_stoich)
            reaction.add_metabolites(charge_stoich)
            transported_protons, transported_charge,outer_h,outer_z,inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_dp,rxn_already_balanced=False)



        #reaction_swapped, swapped, inner_comp, outer_comp = _transport_direction(reaction)

        pH_inner = reaction.model.pH[inner_comp]
        pH_outer = reaction.model.pH[outer_comp]
        e_potential_difference = reaction.model.phi_dict[inner_comp][outer_comp] #ToDo check these are the correct way around? 
        dg_protons = (
            transported_protons
            * R
            * reaction.model.T
            * np.log(10.0)
            * (pH_inner - pH_outer)) #swapped flag ensures inner and outer compartments match original reaction direction 

        dg_electrostatic = FARADAY * transported_charge * e_potential_difference
            
    else:
        dg_protons = Q_(0,'kJ/mol')
        dg_electrostatic = Q_(0,'kJ/mol')
        
    return -dg_protons-dg_electrostatic, -dg_protons, -dg_electrostatic

def leading_zeros(decimal):
    '''Calculate the leading 0s in a decimal i.e 0.001 returns 2'''
    return False if abs(decimal) <= 1e-6 else -floor(log10(abs(decimal))) - 1

def net_elements(reaction, balance_mg = True, round_dp=False, rxn_already_balanced = True):
    """Calculates the net protons, charge, and mg of a reaction"""
    if balance_mg is False:
        pMg = Q_(14,'')
    else:
        pMg = None
    if reaction.boundary == True:
        return {}, {}
    
    net_elements = {}
    # transported_free_h = 0
    #for non transporters 
    if len(reaction.compartments) < 2 and reaction.boundary is False:
        net_protons = 0
        net_charge = 0
        net_mg = 0

        for metabolite, stoich in reaction.metabolites.items():
            charge, protons, mg, microspecies_df = metabolite.average_charge_protons(pMg=pMg, round_dp=round_dp,cobra_formula=reaction._cobra_formula)           
            net_protons += (protons*stoich)
            net_charge += (charge*stoich)
            net_mg += (mg*stoich)

        comp = list(reaction.compartments)[0]
        net_elements[reaction.model.proton_dict[comp]] = net_protons
        net_elements[reaction.model.charge_dict[comp]] = net_charge
        net_elements[reaction.model.mg_dict[comp]] = net_mg

        return net_elements, {}

    #for transport reactions
    elif len(reaction.compartments) == 2:

        n_h_inner, z_inner, n_h_outer, z_outer, inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_dp)

        #catch new round defined 
        if balanced is False:
            round_estimate = leading_zeros(n_h_inner + n_h_outer + 1e-6) #estimate what the previous rounding was based on the remainder
            n_h_inner, z_inner, n_h_outer, z_outer, inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_estimate, rxn_already_balanced=rxn_already_balanced)
            #this does not work if new round is less than old round
            #in this case just 
            i = 0 
            while balanced is False:
                round_estimate +=1
                i+=1
                n_h_inner, z_inner, n_h_outer, z_outer, inner_comp, outer_comp, balanced = transported_c_h(reaction, round_dp = round_estimate, rxn_already_balanced=rxn_already_balanced)

                if i >=10: #try adjusting round 10 times then give up
                    break
                

        reaction_swapped, swapped, inner_comp, outer_comp = _transport_direction(reaction)
        transported_mets = calc_transported_mets(reaction)

        # define compartment conditions
        # for transported mets set low Mg conc as we assume Mg is not transported 
        pH_inner = reaction.model.pH[inner_comp]
        ionic_strength_inner = reaction.model.I[inner_comp]
        temperature = reaction.model.T

        net_protons_inner = 0
        net_protons_outer = 0
        net_mg_inner = 0
        net_mg_outer = 0

        if not balanced:
            free_prot_trans_inner = 0
            free_prot_trans_outer = 0
        else:
            free_prot_trans_inner = n_h_inner
            free_prot_trans_outer = n_h_outer

        free_mg_trans_outer = 0
        free_mg_trans_inner = 0

        net_reac_protons_inner = 0
        net_reac_protons_outer = 0
        net_reac_mg_inner = 0
        net_reac_mg_outer = 0

       # print(free_prot_trans_inner, free_prot_trans_outer, 'free_prot_trans')

        #deal with transported metabolites first
        for met, stoich in transported_mets.items(): # this includes user defined transported metabolites
            average_charge, average_protons, average_mg, ms_df = met.average_charge_protons(pMg = pMg, round_dp=round_dp,cobra_formula=reaction._cobra_formula) 
          #  print(met.id, average_charge, average_protons)
            if average_mg is None:
                average_mg = 0
            #define major microspeceis of transported metabolites assuming inner compartment
            major_charge, major_protons, major_Mg = major_microspecies(met, pH_inner, Q_(14,''), ionic_strength_inner, temperature,cobra_formula=reaction._cobra_formula) #consider inner conditions for transported metabolites 
          #  print(met.id, major_charge, major_protons)

            if met in reaction.model.mg_dict.values():#special case for magnesium transporter - otherwise normally magnesium does not cross membrane with metabolites
                major_mg = 1
                major_protons = 0
                major_charge = 2
            else:
                major_mg = 0

            # if the metabolite is in the user defined transported mets and this is in the inner comp then do as usual
            # but also subtract this from the other compartment as we assume it is transported  
            # if the metabolite is in the user defined transported mets and this is in the outer comp then do as susual 
            # but also subtract from the other compartmen 
           #print(met.id, free_prot_trans_outer, free_prot_trans_inner, 'prot_inner_outer')
            if met.compartment == inner_comp: 
                #calculate protons lost/gained converting from major microspecies to average microspecies and vice versa
                net_protons_inner += ((major_protons - average_protons)*stoich)
                net_mg_inner += ((major_mg- average_mg)*stoich)

                if balanced: 
                    #calculate free protons crossing membrane 
                    if  met not in reaction.model.proton_dict.values():
                        free_prot_trans_inner -= stoich*major_protons #after calculating net protons crossing membrane any that are not part of the transported metabolite must be free protons transported
                        # special case of transport and reaction 
                        # transported metabolites are just deinfed once in one compartment so we also need to account for protons on other compartment
                        if ((met in reaction.transported_mets) if reaction.transported_mets is not None else False):
                           # print(met.id, free_prot_trans_outer,stoich*major_protons ,'free protons crossing')

                            free_prot_trans_outer += stoich*major_protons #after calculating net protons crossing membrane any that are not part of the transported metabolite must be free protons transported
                            #print(met.id, free_prot_trans_outer,stoich*major_protons ,'free protons crossing')
                else:
                    if  met in reaction.model.proton_dict.values():
                        free_prot_trans_inner += stoich 

                #calculate free mg crossing membrane 
                if met in reaction.model.mg_dict.values():
                    free_mg_trans_inner += stoich

            if met.compartment == outer_comp:
              #  print(met.id, outer_comp) 
                #calculate protons lost/gained converting from major microspecies to average microspecies and vice versa
                net_protons_outer += (major_protons - average_protons)*stoich
                net_mg_outer += (major_mg - average_mg)*stoich
             #   print(net_protons_outer, 'net_protons_outer')
                if balanced:
                    #calculate free protons crossing membrane 
                    if  met not in reaction.model.proton_dict.values():
                   #     print(free_prot_trans_outer, 'free_prot_trans_outer',stoich,major_protons)

                        free_prot_trans_outer -= stoich*major_protons #ToDo maybe change this to average protons?
                     #   print(free_prot_trans_outer, 'free_prot_trans_outer')
                        # special case of transport and reaction 
                        # transported metabolites are just deinfed once in one compartment so we also need to account for protons on other compartment
                        if ((met in reaction.transported_mets) if reaction.transported_mets is not None else False):
                          #  print(met.id, free_prot_trans_inner,stoich*major_protons ,'free protons crossing inner')

                            free_prot_trans_inner += stoich*major_protons #after calculating net protons crossing membrane any that are not part of the transported metabolite must be free protons transported
                         #   print(met.id, free_prot_trans_inner,stoich*major_protons ,'free protons crossing inner')

                else:
                #calculate free protons crossing membrane 
                    if met in reaction.model.proton_dict.values():
                        free_prot_trans_outer += stoich

                #calculate free mg crossing membrane 
                if met in reaction.model.mg_dict.values():
                    free_mg_trans_outer += stoich            

          #  print(net_protons_inner, net_protons_outer, 'net_protons' )
          #  print(free_prot_trans_inner, free_prot_trans_outer, 'trans_protons' )

        
        # now deal with chemical reaction part 
        # this is reaction - transported metabolites


        for met, stoich in reaction.metabolites.items():

            #if the metabolite is in the predefined transported and reaction metabolites then treat it as a reaction metabolite
            if ((met in reaction.transported_mets) if reaction.transported_mets is not None else False):
                stoich = stoich 
             #   print(met.id, stoich, 'a')
            #otherwise it might be a normal reaction metabolite
            elif met in transported_mets:
                #for metabolites that are separate reaction and transport metabolites e.g.Pi in PiABC subtract the transport part 
                stoich -= transported_mets[met] 
             #   print(met.id, stoich, 'b')


            charge, protons, mg, ms_df = met.average_charge_protons(pMg = pMg, round_dp=round_dp,cobra_formula=reaction._cobra_formula) 
            if mg == None:
                mg = 0

            # special case of chemical reaction at same time as transport 
            # consider major microspecies - defined by inner compartment 
            if ((met in reaction.transported_mets) if reaction.transported_mets is not None else False):
                #define major microspeceis of transported metabolites assuming inner compartment
                charge, protons, mg = major_microspecies(met, pH_inner, Q_(14,''), ionic_strength_inner, temperature,cobra_formula=reaction._cobra_formula) #consider inner conditions for transported metabolites 

                if met not in reaction.model.proton_dict.values():
                    if met not in reaction.model.mg_dict.values():
                        # note swapping of inner and outer compartments
                        # this assumes a metabolite invovled in simultaneous transport and reaction is transported and then reacts
                        if met.compartment == outer_comp: 
                            net_reac_protons_inner += stoich*protons
                            net_reac_mg_inner += stoich*mg

                        if met.compartment == inner_comp:
                            net_reac_protons_outer += stoich*protons
                            net_reac_mg_outer += stoich*mg
            else:             
                #net protons and mg for non-transport parts of reaction 
                #ignore protons and mg here as we are recalculating the proton balance of the non-transport part of reaction
                if met not in reaction.model.proton_dict.values():
                    if met not in reaction.model.mg_dict.values():

                        if met.compartment == inner_comp:
                            net_reac_protons_inner += stoich*protons
                            net_reac_mg_inner += stoich*mg

                        if met.compartment == outer_comp:
                            net_reac_protons_outer += stoich*protons
                            net_reac_mg_outer += stoich*mg

              #  print(met.id, protons, 'x')

       # print(net_reac_protons_inner,free_prot_trans_inner,net_protons_inner, 'y' )
      #  print(net_reac_protons_outer,free_prot_trans_outer,net_protons_outer, 'z')

        #if the transported_c_h did not return a balanced reaction then we must account for this when calculaitng charge
        if not balanced:
         #   print(net_protons_inner , -net_reac_protons_inner ,-2*net_reac_mg_inner , 2*net_mg_inner, 'a')
         #   print(net_protons_outer , -net_reac_protons_outer ,-2*net_reac_mg_outer , 2*net_mg_outer, 'b')

            z_inner += (net_protons_inner - net_reac_protons_inner -2*net_reac_mg_inner + 2*net_mg_inner)
            z_outer += (net_protons_outer - net_reac_protons_outer -2*net_reac_mg_outer + 2*net_mg_outer)

    
        #if reaction has already defined transported protons and charge then use them here instead of calculated values  
        if reaction._transported_h is not None:

            z_inner -= free_prot_trans_inner #remove transported free proton component from calcualted charge  
            z_outer -= free_prot_trans_outer

            #update with new value 
            free_prot_trans_inner = reaction._transported_h[inner_comp]
            free_prot_trans_outer = reaction._transported_h[outer_comp]

            z_inner += free_prot_trans_inner #add new transported free proton charge  
            z_outer += free_prot_trans_outer      

        transported_free_h = {inner_comp:free_prot_trans_inner,
                        outer_comp:free_prot_trans_outer}

        #ToDo: check if this is already accoutned for in transported_c_h function 
       # if reaction.transported_charge is not None:
        #    free_charge_trans_inner = reaction._transported_charge[inner_comp]
         #   free_charge_trans_outer = reaction._transported_charge[outer_comp]
#
 #           z_inner += free_charge_trans_inner#add new transported free charge  
  #          z_outer += free_charge_trans_outer 


        #calculate total proton and mg balance 
        
        net_inner_protons= -1*net_reac_protons_inner+free_prot_trans_inner+net_protons_inner
        net_outer_protons = -1*net_reac_protons_outer+free_prot_trans_outer+net_protons_outer

        net_inner_mg = -1*net_reac_mg_inner+net_mg_inner+free_mg_trans_inner
        net_outer_mg = -1*net_reac_mg_outer+net_mg_outer+free_mg_trans_outer

        net_elements = {}
        net_elements[reaction.model.proton_dict[inner_comp]] = net_inner_protons*-1 #ToDO check if -1 is needed here?
        net_elements[reaction.model.charge_dict[inner_comp]] = z_inner *-1 
        net_elements[reaction.model.mg_dict[inner_comp]] = net_inner_mg*-1

        net_elements[reaction.model.proton_dict[outer_comp]] = net_outer_protons*-1
        net_elements[reaction.model.charge_dict[outer_comp]] = z_outer *-1
        net_elements[reaction.model.mg_dict[outer_comp]] = net_outer_mg*-1

        

        for key, value in net_elements.items():
            net_elements[key] = round(value, 6) #ToDo is there a better way for dealing with floaitng point issues?

    return net_elements, transported_free_h


def reaction_balance(reaction, balance_charge = True, balance_mg = True, round_dp = False, rxn_already_balanced=False):
    ''' 
    Docstring for reaction_balance

    :param reaction: Description
    :type reaction: 
    :param balance_charge: Description
    :type balance_charge: 
    :param balance_mg: Description
    :type balance_mg: 
    :param round_dp: Description
    :type round_dp: 
    :param rxn_already_balanced: Description
    :type rxn_already_balanced: 
    :return: Description
    :rtype: Any'''

    transporter = False

    if reaction.boundary == False:
        if len(reaction.compartments) == 3:
            list_sub_rxn = reaction.split_reaction() # this returns a list of sub reactions (2compartments) and modifies the input reaction
            reaction.model.charge_dict
            reaction.model.proton_dict
            reaction.model.mg_dict
            if len(reaction.compartments)>0:
                assert len(reaction.compartments) <3, 'The reaction still has more than 2 compartments' # this means that one metabolite is not conserved over the transport of 2 compartments, and we cannot deal we this case yet
                a= reaction_balance(reaction, balance_charge=True, balance_mg=False)#maybe the reaction has been totally emptied
                drg_transport, dg_protons, dg_electrostatic = calc_drGtransport(reaction,  round_dp = round_dp)
            else:
                drg_transport, dg_protons, dg_electrostatic = Q_(0,'kJ/mol'), Q_(0,'kJ/mol'), Q_(0,'kJ/mol')
            for sub_rxn in list_sub_rxn:
                reaction.model.add_reactions([sub_rxn])
                a=reaction_balance(sub_rxn, balance_charge=True, balance_mg=False)
                reaction.add_metabolites(sub_rxn.metabolites)
                sub_drg_transport, sub_dg_protons, sub_dg_electrostatic = calc_drGtransport(sub_rxn,  round_dp = round_dp)
                drg_transport += sub_drg_transport
                dg_protons += sub_dg_protons
                dg_electrostatic += sub_dg_electrostatic
                #tuples don't support augmented assignment 
                sub_rxn.remove_from_model()
            #set the drGtransport of the reaction to the sum of the sub reactions and the core one, because calc_drG_transport is also native for 2 compartments reaction
            reaction.drGtransport = drg_transport
            reaction.drG_h_transport = dg_protons
            reaction.drG_c_transport = dg_electrostatic
            # print(reaction.reaction,'3 compartments reaction balanced and drGtransport added', reaction.drGtransport)
            return
        elif len(reaction.compartments) ==2:
            transporter = True
            #catch transporters with unknown transported metabolites
            if (calc_transported_mets(reaction) == {}) and (reaction.transported_mets) is None and (reaction._transported_h) is None:
                #raise exception 
                warnings.warn(f"WARNING: {reaction.id} is a transporter but the transported metabolites could not be automatically calculated and no transported metabolites were manually defined, please check reaction stoichiometry or define reaction.transported_mets")
                reaction.report_prio=1
        net_elements, transported_free_h = thermo_flux.tools.drg_tools.net_elements(reaction, balance_mg,  round_dp = round_dp)

        #if the reaction is a transporter and transported_free_h do not cancel out then raise a warning
        if transporter:
            if abs(sum(transported_free_h.values()))>1e-6:
                warnings.warn(f"WARNING: {reaction.id} is not balanced and could not be automatically balanced, please check reaction stoichiometry or define reaction.transported_h", stacklevel=2)
                reaction.report_prio=1
                return
        
        #update transported_free_h
        reaction._transported_h = transported_free_h

        #remove all proton charge and mg from reaction before balancing as they will be replaced
        # only for transporters
        if transporter:
            for met, stoich in reaction.metabolites.items():
                if met in reaction.model.proton_dict.values():
                    reaction.add_metabolites({met: -1*stoich})

                if balance_charge: #only balance charges for transport reactions
                    if met in reaction.model.charge_dict.values():
                        reaction.add_metabolites({met: -1*stoich})
                if balance_mg:
                    if met in reaction.model.mg_dict.values():
                        reaction.add_metabolites({met: -1*stoich})


            
        if type(round_dp) is int:
            reaction.add_metabolites({met: -1*round(stoich, round_dp) for met, stoich in net_elements.items() if met in reaction.model.proton_dict.values() })
        else:
            reaction.add_metabolites({met: -1*stoich for met, stoich in net_elements.items() if met in reaction.model.proton_dict.values() })

        if balance_charge:
            if transporter: #only balance charges for transport reactions
                if type(round_dp) is int:
                    reaction.add_metabolites({met: -1*round(stoich, round_dp) for met, stoich in net_elements.items() if met in reaction.model.charge_dict.values() })
                else:
                    reaction.add_metabolites({met: -1*stoich for met, stoich in net_elements.items() if met in reaction.model.charge_dict.values() })

        if balance_mg:
            if type(round_dp) is int:
                reaction.add_metabolites({met: -1*round(stoich, round_dp) for met, stoich in net_elements.items() if met in reaction.model.mg_dict.values() })
            else:
                reaction.add_metabolites({met: -1*stoich for met, stoich in net_elements.items() if met in reaction.model.mg_dict.values() })


        reaction.balanced = True

    return

from math import  copysign
from itertools import product

def generate_combinations(dictionary):
    keys = dictionary.keys()
    values = dictionary.values()
    combinations = product(*values)
    return [dict(zip(keys, combination)) for combination in combinations]


def new_reaction_name(reaction, charge_states):
    '''Function to generate a new reaction name for a transporter variant with a specific charge state.
    
    Parameters
    ----------
    reaction: tmodel.reaction
        reaction to add variants of
    charge_states: list of int
        list of charge states of the transported metabolites
        
    Returns
    -------
        str: name of the new reaction
    '''

    ID_compartment = ''
    if any([reaction.id.endswith('_'+(''.join(reaction.compartments))), 
        reaction.id.endswith('_'+(''.join(reaction.compartments))[::-1])]):    
        ID_compartment = '_'+reaction.id.split("_",100)[-1] #get compartment form metabolite name (splits id at _ with max 100 splits then takes the last string)

    charge_states_string = ''.join(['+' + str(n) if n > 0 else str(n) for n in charge_states])
    name = reaction.id.replace(ID_compartment,"") + charge_states_string + ID_compartment

    return name


def add_transporter_varaints(reaction, add_charge_neutral = True, balance_charge = False, round_dp = False):
    '''Function to add all  transporter variants of a reaction.
    Variants are added for all species of transported metabolites with an abundance of >10% in the inner compartment.
    Additional transporters to represent a charge neutral transporter can also be added (including intermediate charge transport)

    Parameters
    ----------
    reaction: tmodel.reaction
        reaction to add transporter variants of
    add_charge_neutral: bool
        if True then a charge neutral variants of the transporter will be added
    balance_charge: bool
        if True then the reaction will be balanced for charge (if original reaction contains charge then this will be forced True)
    round_dp: int
        for charge and proton balaning if rounded values are used (should match original reaction rounding)

    Returns
    -------
    list of ThermoReaction objects representing the transporter variants
    '''
    
    transporter_variants = set({}) 


    for met in reaction.metabolites:
        if met in reaction.model.charge_dict.values():
            balance_charge = True

    reaction_balance(reaction, balance_charge=balance_charge, balance_mg=False,  round_dp = round_dp)

    n_h_inner, z_inner, n_h_outer, z_outer,inner_comp, outer_comp, balanced = reaction.transported_c_h(round_dp = round_dp)

    #identify transported metabolite

    #identify number of different species of that metabolite 

    #add protons to represent all the different speceis 

    #add all transporter variants i.e one for each possible charge state of the transported metabolite(s) that is >10% abundant
    
    
    #identify the transported metabolites 
    tmets = [met for met in calc_transported_mets(reaction)
            if (met in reaction.model.proton_dict.values())
            & (met not in reaction.model.charge_dict.values())
                & (met.compartment == inner_comp)]

    if len(tmets) == 0:
        tmets = [met for met in calc_transported_mets(reaction)
            if (met in reaction.model.proton_dict.values())
            & (met not in reaction.model.charge_dict.values())]

    #if tmets is empty then there are no transported metabolites other than charge and protons and species variants should not be added
    #this prevents adding variants for things like ATPsynthase or ETC redox reactions
    if len(tmets) == 0:
        return transporter_variants

    #identify the species > 10% in the innter compartment 
    tmets_species = {}
    for met in tmets:
        species = met.average_charge_protons(pMg = Q_(14,''), accuracy=0.1)[3] #accuracy = 0.1 drops species < 10% abundant 
        if species is not None:
            if 'abundance_norm' in species:
                #number of charge states ignoring magneiusm with >10% abundance
                tmets_species[met] = species.loc[(species['number_magnesiums']==0)&(species['abundance_norm']!=0)]['charge'].unique()
            
    #create a new reaction for transport of each possible species of combination thereof
    species_combinations =  generate_combinations(tmets_species)

    new_reactions = []


    #if there's only one set of species combinations then all the relevant species
    # are already represented by assuming just the major species is transported
    #otherwise there are more combinations that should be added

    if len(species_combinations) >1:
        #create new reaction for each possible species >10% and their combinations 
        for species in species_combinations:
            #create a new reaction with the relevant name
            #naming convention is reaction_name, charge_state1, charge_state2, etc.
            new_reaction = reaction.copy()
            name = new_reaction_name(reaction, list(species.values()))

            inner_proton = reaction.model.proton_dict[inner_comp]
            outer_proton = reaction.model.proton_dict[outer_comp]

            #calculate how many protosn to add to represent the required charge states
            # of the transported metabolites

            for met, charge_state in species.items():
                protons_to_add = charge_state - met.major_microspecies(pMg = Q_(14,''))[0]

                met_stoich = new_reaction.metabolites[met]

                new_reaction.add_metabolites({inner_proton: met_stoich*((protons_to_add)), 
                                        outer_proton: -met_stoich*((protons_to_add))})
        
                if new_reaction._transported_h is None:
                    new_reaction._transported_h = {inner_comp:met_stoich*(protons_to_add),
                                                    outer_comp:-met_stoich*(protons_to_add)}
                else:
                    new_reaction._transported_h[inner_comp] += met_stoich*(protons_to_add)
                    new_reaction._transported_h[outer_comp] += -met_stoich*(protons_to_add)


            new_reactions.append(new_reaction)


    #create new reaction for each possible variant up to charge neutral 
    #if z_inner is already 0 then this will add nothing as a charge neutral variant already exists
    if add_charge_neutral is True: 
        for i in range(0,int(abs(z_inner))):
            new_reaction = reaction.copy()

            name = new_reaction_name(reaction, [i])
        
            new_reaction.id = name
            new_reaction.name = name

            inner_proton = reaction.model.proton_dict[inner_comp]
            outer_proton = reaction.model.proton_dict[outer_comp]

            new_reaction.add_metabolites({inner_proton: -1*(copysign(i+1,z_inner)), 
                                        outer_proton: -1*(copysign(i+1,z_outer))})
            
            if new_reaction._transported_h is None:
                new_reaction._transported_h = {inner_comp:-1*copysign(i+1,z_inner),
                                                outer_comp:-1*copysign(i+1,z_outer)}
            else:
                new_reaction._transported_h[inner_comp] += -1*copysign(i+1,z_inner)
                new_reaction._transported_h[outer_comp] += -1*copysign(i+1,z_outer)
            
            new_reactions.append(new_reaction)

        ##add the new reactions to the model


    rxn_in_model = False
    model_rxn_ids = [rxn.id for rxn in reaction.model.reactions]
    
    model = reaction.model ###if we delete reaction from model, then we can't access reaction.model 

    for new_reaction in new_reactions:

        #balance the new reaction for comoparison with the model (model must be balanced for this to work)
        #if the reaction is not in the model but thermo-flux has already added a reaction with this name then delete it
        #note this may delete a reaction already in the model with the same name but different stoichiometry... 

        if new_reaction.id in model_rxn_ids:
            model.reactions.get_by_id(new_reaction.id).delete()

        model.add_reactions([new_reaction])
        reaction_balance(new_reaction, balance_charge=balance_charge, balance_mg=False,  round_dp = round_dp)
        model.reactions.get_by_id(new_reaction.id).delete()

        #if an identical reaction is already in the model then do not replace it 
        for rxn in model.reactions:
            if new_reaction.metabolites == rxn.metabolites:
                print('transporter variant',new_reaction.id, 'already in model as', rxn)
                transporter_variants.add(rxn)
                rxn_in_model = True
                #delete the new reaction we just added as it is already in the model (maybe under a different name)
                break

        if rxn_in_model is False:
            transporter_variants.add(new_reaction)
            model.add_reactions([new_reaction])


    return transporter_variants



def add_transporter_charge_varaint(reaction, charge_state, round_dp = False):
    '''Function to add a transporter variant with a specific net charge transport. Protons are added to either side of
    the equation to represent different charge states. Negative values are represented by protons with an opposite
    transport stoichiometry. 
    
    Parameters
    ----------
    reaction: tmodel.reaction
        reaction to add variants of
    charge_state: int
        desired net charge to be transported (can represent different charge states of the transported metabolite)
    round_dp: int
        for charge and proton balaning if rounded values are used 

    Returns
    -------
    ThermoReaction with a the desired net charge transport

    '''

    
    #firstly ensure the reaciton is balanced 
    reaction_balance(reaction, balance_charge=False, balance_mg=False,  round_dp = round_dp)

    #calcualte the transported charge and protons as the reaction is already written 
    n_h_inner, z_inner, n_h_outer, z_outer,inner_comp, outer_comp, balanced = reaction.transported_c_h(round_dp = round_dp)

    #calculate how many protons to add based on desired charge state
    protons_to_add = charge_state - z_inner

    #if the charge state requested is already the original reaction then do nothing and return the original reaction
    if protons_to_add == 0:
        return reaction

    #how to relate sign of charge change with direction of reaction? 

    inner_proton = reaction.model.proton_dict[inner_comp]
    outer_proton = reaction.model.proton_dict[outer_comp]

    
    #create new reaction 
    new_reaction = reaction.copy()

    ID_compartment = ''
    if any([reaction.id.endswith('_'+(''.join(reaction.compartments))), 
        reaction.id.endswith('_'+(''.join(reaction.compartments))[::-1])]):    
        ID_compartment = '_'+reaction.id.split("_",100)[-1] #get compartment form metabolite name (splits id at _ with max 100 splits then takes the last string)


    name = reaction.id.strip(ID_compartment) + str(int(charge_state)) + ID_compartment
    new_reaction.id = name
    new_reaction.name = name

    #add the new reaction to the model
    #if the reaction is already in the model then delete it
    if new_reaction.id in [rxn.id for rxn in reaction.model.reactions]:
        reaction.model.reactions.get_by_id(new_reaction.id).delete()

    reaction.model.add_reactions([new_reaction])

    new_reaction.add_metabolites({inner_proton: 1*((protons_to_add)), 
                                    outer_proton: -1*((protons_to_add))})
    
    if new_reaction._transported_h is None:
        new_reaction._transported_h = {inner_comp:1*(protons_to_add),
                                        outer_comp:-1*(protons_to_add)}
    else:
        new_reaction._transported_h[inner_comp] += 1*(protons_to_add)
        new_reaction._transported_h[outer_comp] += -1*(protons_to_add)
    
    reaction_balance(new_reaction, balance_charge=False, balance_mg=False,  round_dp = round_dp)

    return new_reaction
           


