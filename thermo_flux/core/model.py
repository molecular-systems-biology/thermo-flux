"""A model with thermodynamic properties"""

from cobra.core import Model, Gene, Group, Reaction, Metabolite
from cobra.util.array import create_stoichiometric_matrix
from thermo_flux.core.metabolite import ThermoMetabolite
from thermo_flux.core.reaction import ThermoReaction

from ..tools import drg_tools
from ..io.output import generate_report

from copy import copy, deepcopy

from typing import Dict

from cobra.core.dictlist import DictList

from equilibrator_api import ComponentContribution, Q_

from pathlib import Path

import numpy as np
import pandas as pd

import gurobipy as gp
from ..solver import gurobi

from ..utils.vis import progressbar
import thermo_flux
import warnings
def custom_formatwarning(message, category, filename, lineno, line=None):
    return f'{message}\n'
warnings.formatwarning = custom_formatwarning
from warnings import warn

default_rmse_inf = Q_(1e5, 'kJ/mol')
default_max_drG = Q_(1e8, 'kJ/mol')
default_T = Q_(303.15, 'kelvin')
default_pH = Q_(7)
default_phi = Q_(0, 'volt')
default_I = Q_(0.2,'M')
default_pMg = Q_(3)


class ThermoModel(Model):
    """Thermodynamic model object.

    takes cobra model and adds thermodynamic properties

    """

    def __init__(
        self,
        model,
        pH: Dict[str, Q_] = None,
        I: Dict[str, Q_] = None,
        pMg: Dict[str, Q_] = None,
        T: Dict[str, Q_] = None,
        phi: Dict[str, Q_] = None,
        gdiss_lim: Q_ = None,
        dfG0_cov_sqrt=None,
        update_thermo_info=False,
        rmse_inf=default_rmse_inf,
        cc=None,
        lc=None,
        m=None,
        max_drG=None,
        split_biomass=True,
        update_biomass_dfG0=False,
        add_charge_exchange=True,
        phi_dict :Dict[str,Q_] = None,
        **kwargs):

        Model.__init__(self)
        self._pH = pH
        self._I = I
        self._pMg = pMg
        self._T = T
        self._phi = phi
        self._gdiss_lim = gdiss_lim
        self._proton_dict = None
        self._charge_dict = None
        self._mg_dict = None
        self._rmse_inf = rmse_inf
        self._cc = cc
        self._lc = lc
        self._m = m 
        self._max_drG = max_drG
        self._phi_dict = phi_dict
        self._biomass_rxn=False
        
        for attr, value in model.__dict__.items():
            if attr not in ['metabolites', 'genes', 'reactions']:
                self.__dict__[attr] = value

        self.metabolites = DictList()

        for met in model.metabolites:
            tmet = ThermoMetabolite(met, model=self)
            self.metabolites.append(tmet)

        self.reactions = DictList()
        for reaction in model.reactions:
            treaction = ThermoReaction(reaction, model=self)
            self.reactions.append(treaction)

        if hasattr(model, '_parameters'):
            self.update_parameters()

        if update_thermo_info:
            self.update_thermo_info()

        if split_biomass:
            self.split_biomass()

        if update_biomass_dfG0:
            #biomass_dfG0 calcuation requires split biomass
            if split_biomass is False:
                self.split_biomass()
            self.update_biomass_dfG0()

        if add_charge_exchange:
            self.add_charge_exchange()
           

    def get_compounds(self, search = False, update_annotations = False, check_consisitency = True):
        searched_mets = []
        for met in progressbar(self.metabolites,'', 40, item_label_attribute = 'id'):
            if met._compound is None:
                cpd, annotation, formula, inchi , searched = drg_tools.get_suitable_ids(met,search=search, update_annotations = update_annotations)
                met.annotation.update(annotation)
                met.compound = cpd
                if searched:
                    searched_mets.append(met)

            if check_consisitency:
                if met.check_consistency() is False:
                    eq_atom_bag = {}
                    if met.compound is not None:
                        if met.compound.atom_bag is not None:
                            eq_atom_bag = met.compound.atom_bag
                            warnings.warn('WARNING: {} has inconsistent formula with eQuilibrator identified metabolite. COBRA:{}, eQuilibrator:{}. Defining as unknown compound'.format(met.id, met.elements, eq_atom_bag), stacklevel=2)
                            met.unknown = True

        return searched_mets


    def update_thermo_info(self, fit_unknown_dfG0=False, search = False, round_dp = False, report=False):

        print('Identifying compounds...')
        self.get_compounds(search = search)        
        #load charge, proton and mg dicts in case these metabolites are not present
        self.charge_dict
        self.proton_dict
        self.mg_dict

        print("Estimating dfG0'...")
        dfG0_mean, dfG0prime_mean, dfG0_cov_sqrt, unknowns_basis = drg_tools.calc_dfG0prime(self, fit_unknown_dfG0=fit_unknown_dfG0)
        dfG_SE = [SE for SE in list(np.sqrt((dfG0_cov_sqrt @ dfG0_cov_sqrt.T).diagonal()).flat)]
        
        for i, met in enumerate(self.metabolites):
            met.dfG0 = dfG0_mean[i]
            met.dfG0prime = dfG0prime_mean[i]
            met.dfG_SE = dfG_SE[i]
            cpd=met.compound

            #if (cpd is None) or (cpd.can_be_transformed() is False) or (cpd.atom_bag is None) or (cpd.atom_bag == {}):
             #   warn(f"Warning: Metabolite {met.id} doesn't have an eQuilibrator formula.")
        
        print("Estimating drG0'...")
        drG0_mean = drg_tools.calc_model_drG0(self)
        drG0prime_mean = drg_tools.calc_model_drG0prime(self)
        drGmprime_mean = drG0prime_mean + drg_tools.calc_phys_correction(self)

        S = create_stoichiometric_matrix(self)
        standard_dgr_Q = (S.T @ dfG0_cov_sqrt)  # sqrt of covariance 
        
        drG_SE = [SE for SE in list(np.sqrt((standard_dgr_Q @ standard_dgr_Q.T).diagonal()).flat)]

        for i, rxn in enumerate(progressbar(self.reactions,'', 40, item_label_attribute = 'id')):
            rxn.drG0 = drG0_mean[i]
            rxn.drG0prime = drG0prime_mean[i]
            rxn.drG_SE = drG_SE[i]
            rxn.drGmprime = drGmprime_mean[i]
        

            if len(rxn.compartments) > 1:

                #if the reaction is a transporter but no transported metabolites are identified then raise a warning
                if thermo_flux.tools.drg_tools.calc_transported_mets(rxn) == {}:
                    warnings.warn(f"WARNING: {rxn.id} is a transporter but no transported metabolites could be automatically identified. Please update the reaction.transported_mets", stacklevel=2)
                    rxn.report_prio=1
                drg_transport, dg_protons, dg_electrostatic = drg_tools.calc_drGtransport(rxn,  round_dp = round_dp)
                rxn.drGtransport = drg_transport
                rxn.drG_h_transport = dg_protons
                rxn.drG_c_transport = dg_electrostatic

            else:
                rxn.drGtransport = Q_(0, 'kJ/mol')
                rxn.drG_h_transport = Q_(0, 'kJ/mol')
                rxn.drG_c_transport = Q_(0, 'kJ/mol')

        self._dfG0_cov_sqrt = dfG0_cov_sqrt
        standard_dgr_Q=standard_dgr_Q.m
        #drop small coefficents to help with numerics
        standard_dgr_Q[np.abs(standard_dgr_Q) < 1e-5] = 0 

        #drop degrees of freedom that are all 0 
        standard_dgr_Q = standard_dgr_Q[:,np.where(standard_dgr_Q.any(axis=0))[0]]

        self._drG0_cov_sqrt = standard_dgr_Q

        if report ==True:
            return generate_report(self)
        

    def update_parameters(self):

        if self._parameters is not None:
            for key, value in self._parameters.items():
                if 'pH' in key:
                    ph_dict = {}
                    for comp_names, pH in value.items():
                        comp = comp_names.split('[')[1][0]
                        ph_dict[comp] = Q_(pH, '')
                    self._pH = ph_dict
                    
                if 'Ionic Strength' in key:
                    I_dict = {}
                    for comp_names, I in value.items():
                        comp = comp_names.split('[')[1][0]
                        I_dict[comp] = Q_(I, 'M')
                    self._I = I_dict
                    
                if 'Temperature' in key:
                    self._T = Q_(value, 'K')
                    
                if 'lectrical membrane potential' in key:
                    phi_dict = {}
                    for comp_names, phi in value.items():
                        comp = (comp_names.split('|')[0] +
                                comp_names.split('|')[2]).strip(' ')
                        phi_dict[comp] = Q_(phi, 'V')
                    self._phi = phi_dict

    # def biomass(self):
    #     '''identify a biomass reaction'''
    #     # identify any existing biomass reactions
    #     for rxn in self.reactions:
    #         if 'biomas' in rxn.id.lower(): #only set it once
    #             self._biomass_rxn = rxn
    #     self._biomass_rxn.remove_from_model()
                


    def split_biomass(self, biomass_rxn = None):
        '''Automatically expand the biomass reaction to include transport of
            biomass out of the cell. This ensures the gibbs energy balance is
            maintained as biomass contains protons that have been transported
            into the cell, the energy of those protons leaving the cell must
            also be accounted for.'''

        

        # create a biomass metabolite
        biomass_c = ThermoMetabolite(Metabolite('biomass_c',
                                                name="biomass",
                                                compartment='c',
                                                charge=0))
        
        biomass_e = ThermoMetabolite(Metabolite('biomass_e',
                                                name="biomass",
                                                compartment='e',
                                                charge=0))
        
        biomass_c.biomass = True
        biomass_e.biomass = True

        # add biomass metabolite to biomass reaction
        if self._biomass_rxn == False:
            # self.reactions.biomass_c.remove_from_model()
            for rxn in self.reactions:
                if 'biomas' in rxn.id.lower(): #only set it once ?
                    # rxn.id = 'biomass_c'
                    rxn._ignore_snd = True
                    rxn.add_metabolites({biomass_c: 1.0})
            self._biomass_rxn=True
        self.repair()

        # create a biomass transport reaction
        reaction = ThermoReaction(Reaction(id="biomass_ce"))
        
        reaction.lower_bound = -1*reaction.upper_bound

        reaction.add_metabolites({
            biomass_c: -1.0,
            biomass_e: 1.0,
            })
        reaction._ignore_snd = True

        # create an exchange reaction for biomass
        biomassEX = ThermoReaction(Reaction(id="biomass_EX"))
        biomassEX.lower_bound = -1*biomassEX.upper_bound

        biomassEX.add_metabolites({biomass_e: -1.0})
        biomassEX._ignore_snd = True

        self.add_reactions([reaction, biomassEX])
        

        print('added reaction: ', self.reactions.get_by_id('biomass_ce'))
        print('added reaction: ', self.reactions.get_by_id('biomass_EX'))

    def get_proton_transporters(self):
        """
        Get the proton transporters in the model
        """
        protons = self.proton_dict.values()
        charges = self.charge_dict.values()
        proton_transporters = []
        for rxn in self.reactions:
            rxn_met_list = [met in protons for met in rxn.metabolites if met not in charges]
            if len(rxn_met_list) > 0:
                if all(rxn_met_list):
                    proton_transporters.append(rxn)
        return proton_transporters

    def get_charge_transporters(self):
        """
        Get the charge transporters in the model
        """
        charges = self.charge_dict.values()

        charge_transporters = []
        for rxn in self.reactions:
            if all([met in charges
                    for met in rxn.metabolites]):
                charge_transporters.append(rxn)

        return charge_transporters

    def add_charge_exchange(self):
        '''function to add charge transport and exchange
          reactions to the model. Skips compartment if there is already a proton transport reaction'''
        
        # compartmetns to skip as reaction already exists
        comp_skip = []
        
        # identify if charge transport or exchange reactions already exist
        for rxn in self.reactions:
            if all([met in rxn.model.charge_dict.values()
                    for met in rxn.metabolites]):

                comp_skip.append([comp for comp
                                  in rxn.compartments
                                  if comp != 'c'][0])

        for comp, charge in self.charge_dict.items():
            # for each charge add a transport reaction to the cytosol
            if comp != 'c':
                if comp not in comp_skip:
                    
                    rxn_id = 'charge_c'+comp

                    charge_transport_rxn = ThermoReaction(Reaction(rxn_id))

                    charge_transport_rxn.add_metabolites({charge:-1,
                                                        self.charge_dict['c']:1})
                    
                    charge_transport_rxn._transported_charge = {comp: -1, 'c': 1}

                    charge_transport_rxn.lower_bound = -1*charge_transport_rxn.upper_bound

                    self.add_reactions([charge_transport_rxn])
                    
                    print('added reaction: ',
                           self.reactions.get_by_id(rxn_id))

                    if comp == 'e':
                        charge_transport_rxn = ThermoReaction(Reaction('EX_charge'))

                        charge_transport_rxn.add_metabolites({charge:-1})
                        charge_transport_rxn._transported_charge = {comp: -1}

                        charge_transport_rxn.lower_bound = -1*charge_transport_rxn.upper_bound

                        self.add_reactions([charge_transport_rxn])

                        print('added reaction: ',
                               self.reactions.get_by_id('EX_charge'))

    def add_proton_exchange(self):
        '''function to add charge transport and exchange
          reactions to the model. Skips compartment if there is already a proton transport reaction'''
        
        # compartmetns to skip as reaction already exists
        comp_skip = []
        
        # identify if charge transport or exchange reactions already exist
        for rxn in self.reactions:
            if all([met in rxn.model.proton_dict.values()
                    for met in rxn.metabolites]):

                comp_skip.append([comp for comp
                                  in rxn.compartments
                                  if comp != 'c'][0])

        for comp, proton in self.proton_dict.items():
            # for each charge add a transport reaction to the cytosol
            if comp != 'c':
                if comp not in comp_skip:
                    
                    rxn_id = 'H_c'+comp

                    proton_transport_rxn = ThermoReaction(Reaction(rxn_id))

                    proton_transport_rxn.add_metabolites({proton:-1,
                                                           self.proton_dict['c']:1})
                    
                    proton_transport_rxn._transported_h = {comp: -1, 'c': 1}

                    proton_transport_rxn.lower_bound = -1*proton_transport_rxn.upper_bound

                    self.add_reactions([proton_transport_rxn])
                    
                    print('added reaction: ',
                           self.reactions.get_by_id(rxn_id))

                    if comp == 'e':
                        proton_transport_rxn = ThermoReaction(Reaction('EX_H'))

                        proton_transport_rxn.add_metabolites({proton:-1})
                        proton_transport_rxn._transported_h = {comp: -1}

                        proton_transport_rxn.lower_bound = -1*proton_transport_rxn.upper_bound

                        self.add_reactions([proton_transport_rxn])

                        print('added reaction: ',
                               self.reactions.get_by_id('EX_H'))
                        
    def update_biomass_dfG0(self):
        '''update the dfG0 of the biomass reaction based on the formula of biomass_c
            if no formula is defined it will be automatically calcualted from the biomass_c reaction'''

        #check biomass formula for biomass_c and biomass_e match 
        if self.metabolites.biomass_c.formula != self.metabolites.biomass_e.formula:
            warn('biomass_c and biomass_e formulas do not match using biomass_c formula')
            self.metabolites.biomass_e.formula = self.metabolites.biomass_c.formula

        #calculate biomass formation energy based on the formula of biomass_c
        dfG0_bm = drg_tools.calculate_biomass_dfG0(self.metabolites.biomass_c)

        #if formula has been recaluclated make sure biomass_e has same formula as biomass_c
        self.metabolites.biomass_e.formula = self.metabolites.biomass_c.formula

        #update dfG0 of biomass
        #dfG0prime is calcualted when update_thermo_info is called based on the H in the formula
        self.metabolites.biomass_c.dfG0 = dfG0_bm
        self.metabolites.biomass_e.dfG0 = dfG0_bm

    def solution(self):
        '''return a thermo optimization solution object for the model'''
        return gurobi.get_solution(self)


    @property
    def lc(self):
        if (self._lc is None):
            from equilibrator_assets.local_compound_cache import LocalCompoundCache
            lc = LocalCompoundCache()

            cache_name = self.id +'_compound.sqlite'

            #try and load a local cache with the model id
            if Path(cache_name).is_file(): 
                lc.load_cache(cache_name)
                self._lc = lc

            else: #otherwise generate a new local cache 
                lc.generate_local_cache_from_default_zenodo(cache_name, force_write=True)
                lc.load_cache(cache_name)

                self._lc = lc

        return self._lc

    @lc.setter
    def lc(self, value):
        '''load a local compound cache either as an lc object or from a file location string'''
        from equilibrator_assets.local_compound_cache import LocalCompoundCache
        if type(value) is LocalCompoundCache:
            self._lc = value
        elif type(value) is str:
            lc = LocalCompoundCache()
            lc.load_cache(value)

            self._lc = lc

    @property
    def cc(self):
        if (self._cc is None):
            print('Initializing component contribution object...')
            self._cc = ComponentContribution(rmse_inf=self.rmse_inf, ccache=self.lc.ccache)
        return self._cc

    @cc.setter
    def cc(self, value):
        self._cc = value

    @property
    def rmse_inf(self):
        return self._rmse_inf

    @rmse_inf.setter
    def rmse_inf(self, value):
        if value != self.rmse_inf:
            self._rmse_inf = value
            self.cc.predictor.preprocess.RMSE_inf = value.m_as('kJ/mol')
            self.update_thermo_info()

    @property
    def pH(self):
        if self._pH is None:
            pH_dict = {}
            for comp in self.compartments:
                pH_dict[comp] = default_pH
            self._pH = pH_dict

        if self._pH.keys() != self.compartments.keys():
            for comp in self.compartments.keys():
                if comp not in self._pH:
                    self._pH[comp] = default_pH

        return self._pH

    @pH.setter
    def pH(self, value):
        self._pH = value

    @property
    def I(self):
        if self._I is None:
            I_dict = {}
            for comp in self.compartments:
                I_dict[comp] = default_I
            self._I = I_dict

        if self._I.keys() != self.compartments.keys():
            for comp in self.compartments.keys():
                if comp not in self._I:
                    self._I[comp] = default_I
        return self._I

    @I.setter
    def I(self, value):
        self._I = value

    @property
    def pMg(self):
        if self._pMg is None:
            pMg_dict = {}
            for comp in self.compartments:
                pMg_dict[comp] = default_pMg
            self._pMg = pMg_dict

        if self._pMg.keys() != self.compartments.keys():
            for comp in self.compartments.keys():
                if comp not in self._pMg:
                    self._pMg[comp] = default_pMg
        return self._pMg

    @pMg.setter
    def pMg(self, value):
        self._pMg = value

    @property
    def T(self):
        if self._T is None:
            self._T = default_T
        return self._T

    @T.setter
    def T(self, value):
        self._T = value

    @property
    def phi(self):
        if self._phi is None:
            self._phi = {}
            for comp in self.compartments.keys():
                for comp2 in self.compartments.keys():
                    if comp!=comp2:
                        self._phi[comp+comp2] = default_phi
        
        return self._phi

    @phi.setter
    def phi(self, value):
        for key in value.keys():
            if key[::-1] in value:
                if value[key[::-1]] != -1*value[key]:
                    raise ValueError("Duplicate membrane potentials defined that do not match: ({}, {}).".format((key,value[key]), (key[::-1],value[key[::-1]])))

        self._phi = value
        #reset the phi_dict so it is re-populated with any new values
        self._phi_dict = None

        
    @property
    def phi_dict(self):
        if self._phi_dict is None:
            for key in self.phi.keys():
                if key[::-1] in self._phi:
                    if self._phi[key[::-1]] != -1*self._phi[key]:
                        raise ValueError("Duplicate membrane potentials defined that do not match: ({}, {}).".format((key,self._phi[key]), (key[::-1],self._phi[key[::-1]])))

            membrane_pot_rev = {}
            membrane_pot = self.phi.copy()
            for pair, value in membrane_pot.items():
                if len(pair)==2:
                    membrane_pot_rev[(pair[::-1])] = value*-1
                elif len(pair)==4:
                    membrane_pot_rev[(pair[2:]+pair[0:2])] = value*-1
                else:#this means one compartment is composed of 2 letters
                    #take the 2 last letters of the pair and place them in front of the first letter
                    membrane_pot_rev[(pair[1:]+pair[0])] = value*-1

            membrane_pot.update(membrane_pot_rev)
            #Empty del_psi dict with all compartments
            del_psi_dict = {comp : {comp:Q_(0,'V')
                                    for comp in self.compartments.keys()}
            for comp in self.compartments.keys()}
            # print(del_psi_dict)
            to_pop = []
            for k, v in membrane_pot.items():
                if len(k) ==2:
                    # Each k is a two letter string, 1st letter - comp#1, 2nd letter - comp #2.
                    del_psi_dict[k[0]][k[1]]= v
                elif len(k) == 4:
                    # each k is a four letter string, 2 first letter - comp#1, 2last letter - comp #2
                    del_psi_dict[k[0:2]][k[2:4]]= v
                else:
                    #k is a three letter string, we have to find two compartments out of them 
                    if (k[0] in self.compartments and k[1:] in self.compartments and k[0:2] in self.compartments and k[-1] in self.compartments):
                        raise Warning('cannot differentiate membrane potential between {},{} and {},{}'.format(k[0],k[1:],k[0:2],k[-1]))
                    c1= k[0] if k[0] in self.compartments and k[1:] in self.compartments else k[0:2]
                    c2= k[1:] if k[0] in self.compartments and k[1:] in self.compartments else k[-1]
                    if c1 in self.compartments and c2 in self.compartments:
                        del_psi_dict[c1][c2]= v
                    else:
                        to_pop.append(k)
            if len(to_pop)>0:#remove the keys that are not compartments, can't do that during an iteration over the dict
                for k in to_pop:
                    membrane_pot.pop(k)
                    # print('removed {} from membrane_pot dict'.format(k))
            self._phi_dict = del_psi_dict

            return del_psi_dict
        else :
            return self._phi_dict

    @property
    def gdiss_lim(self):
        return self._gdiss_lim

    @gdiss_lim.setter
    def gdiss_lim(self, value):
        self._gdiss_lim = value

    @property
    def dfG0_cov_sqrt(self):
        """Gets the square root of the covariance of the standard formation energies"""
        return self._dfG0_cov_sqrt

    @property
    def proton_dict(self):
        if (self._proton_dict is None) or (self._proton_dict.keys() != self.compartments.keys()):
            proton_dict = drg_tools.proton_dict(self)
            self._proton_dict = proton_dict
        
        return self._proton_dict

    @proton_dict.setter
    def proton_dict(self, value):
        self._proton_dict = value

    @property
    def charge_dict(self):
        if (self._charge_dict is None) or (self._charge_dict.keys() != self.compartments.keys()):
            charge_dict = drg_tools.charge_dict(self)
            self._charge_dict = charge_dict
      
        return self._charge_dict
            
    @charge_dict.setter
    def charge_dict(self, value):
        self._charge_dict = value

    @property
    def mg_dict(self):
        if (self._mg_dict is None) or (self._mg_dict.keys() != self.compartments.keys()):
            mg_dict = drg_tools.mg_dict(self)
            self._mg_dict = mg_dict
        
        return self._mg_dict
            
    @mg_dict.setter
    def mg_dict(self, value):
        self._mg_dict = value
    
    @property
    def m(self):
        if self._m is None:
            self._m = gp.Model("minlp")
        return self._m 

    @m.setter
    def m(self, value):
        if type(value) is gp.Model:
            self._m = value
        else:
            self._m = None

    @property
    def max_drG(self):
        if self._max_drG is None:
            self._max_drG = default_max_drG
        return self._max_drG

    @max_drG.setter
    def max_drG(self, value):      
        self._max_drG = value

    def add_TFBA_variables(self, m = None, conds=[''], error_type='covariance',
                           qnorm=2, alpha=0.95, epsilon=0.5, nullspace=None,
                           gdiss_constraint=False, sigmac_limit=12.3, split_v=False, big_M = False):

        if m is None:
            m = self.m

        m, mvars = gurobi.add_TFBA_variables(self, m=m, conds=conds, error_type=error_type,
                                             qnorm=qnorm, alpha=alpha, epsilon=epsilon, nullspace=nullspace,
                                             gdiss_constraint=gdiss_constraint, sigmac_limit=sigmac_limit,
                                            split_v=split_v, big_M=big_M)

        self.m = m
        self.mvars = mvars


    def regression(self, conds, flux_data, metabolite_data, volume_data,
                   conc_fit=True, flux_fit=True, drG_fit=True, resnorm=1, qm_resnorm = 2,
                   error_type='covariance', conc_units = None,extracellular=None):


        m = self.m
        mvars = self.mvars

        m, mvars = gurobi.regression(self, m=m, mvars= mvars,conds=conds, flux_data=flux_data,
                                     metabolite_data = metabolite_data, volume_data = volume_data,
                                     conc_fit=conc_fit, flux_fit=flux_fit, drG_fit=drG_fit, resnorm=resnorm,
                                    qm_resnorm = qm_resnorm, error_type = error_type, conc_units = conc_units,extracellular=extracellular)

        self.m = m
        self.mvars = mvars


def copy(self) -> "ThermoModel":
        """Modified from COBRApy: Provide a partial 'deepcopy' of the Model.

        All the Metabolite, Gene, and Reaction objects are created anew but
        in a faster fashion than deepcopy.

        Returns
        -------
        cobra.Model: new model copy
        """
        new = self.__class__(self)
        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "genes",
            "notes",
            "annotation",
            "groups"            
        }
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new.__dict__[attr] = self.__dict__[attr]
        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)
        new._pH = deepcopy(self.pH)
        new._I = deepcopy(self.I)
        new._pMg = deepcopy(self.pMg)
        new._T = deepcopy(self.T)
        new._phi = deepcopy(self.phi)
        new._gdiss_lim = deepcopy(self.gdiss_lim)
        new._rmse_inf = deepcopy(self.rmse_inf)
        new._cc = self.cc
        new._max_drG = deepcopy(self.max_drG)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__(metabolite)
            for attr, value in metabolite.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = copy(value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__(reaction)
            for attr, value in reaction.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in reaction._metabolites.items():
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            new_reaction.update_genes_from_gpr()

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group: Group = group.__class__(group.id)
            for attr, value in group.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new
            new.groups.append(new_group)
        for group in self.groups:
            new_group = new.groups.get_by_id(group.id)
            # update awareness, as in the reaction copies
            new_objects = []
            for member in group.members:
                if isinstance(member, ThermoMetabolite):
                    new_object = new.metabolites.get_by_id(member.id)
                elif isinstance(member, ThermoReaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, Reaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, Gene):
                    new_object = new.genes.get_by_id(member.id)
                elif isinstance(member, Group):
                    new_object = new.groups.get_by_id(member.id)
                else:
                    raise TypeError(
                        f"The group member {member!r} is unexpectedly not a "
                        f"metabolite, reaction, gene, nor another group."
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        # it doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new._contexts = []

        return new