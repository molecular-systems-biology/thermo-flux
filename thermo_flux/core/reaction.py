"""A reaction with thermodynamic properties"""

from cobra.core import Reaction
from thermo_flux.core.metabolite import ThermoMetabolite
from equilibrator_api import Q_
import warnings
from thermo_flux.tools import drg_tools 
import pandas as pd

water_annotations = {'kegg': 'C00001',
 'bigg.metabolite': 'h2o',
 'chebi': 'CHEBI:15377',
 'envipath': '5882df9c-dae1-4d80-a40e-db4724271456/compound/969d0227-3069-4e44-9525-7ae7bad84170',
 'hmdb': 'HMDB01039',
 'metacyc.compound': 'OH',
 'reactome': '109276',
 'sabiork.compound': '40',
 'seed': 'cpd00001',
 'metanetx.chemical': 'MNXM114710',
 'synonyms': 'H2O'}


class ThermoReaction(Reaction):
    """Thermodynamic metabolite object.

    takes cobra reaction and adds thermodynamic properties

    """

    def __init__(self,
                 reaction,
                 model=None,
                 drG0:float=0,
                 ignore_snd:bool = False):

        Reaction.__init__(self)

        self._drG0 = drG0
        self._drG0prime = None
        self._drGmprime = None
        self._drGtransport = None
        self._drGtransform = None
        self._drG_h_transport = None
        self._drG_c_transport = None
        self._drG = None
        self._drG_SE = None
        self._transported_h = None
        self._transported_charge = None
        self._transported_mets = None
        self._balanced = False
        self._ignore_snd = ignore_snd
        self._cobra_formula=False

        self._model = model

        for attr, value in reaction.__dict__.items():
            if attr not in ['_model', '_metabolites','_genes']:
                self.__dict__[attr] = value

        self._metabolites = {}
        for met, stoich in reaction._metabolites.items():
            if model is None:
                new_met = ThermoMetabolite(met)
            else:
                new_met = self.model.metabolites.get_by_id(met.id)
            self._metabolites[new_met] = stoich
            new_met._reaction.add(self)

        # ignore 2nd law for boundary reactions
        if self.boundary is True:
            self._ignore_snd = True

        # find specific cases to ignore 2nd law - only if metabolites identified
        # ignore 2nd law for water transport 
        # find reactions where all metabolites are water and they are transporters 
        water_present = [any([water_annotations[identifier] in key if isinstance(key, list) else key == water_annotations[identifier] 
                      for identifier, key in met.annotation.items()
                      if identifier in water_annotations])
                 for met in self.metabolites]

            

        if all(water_present) if bool(water_present) else False: #if reaction is empty then return false and dont set ignore_snd
            self._ignore_snd = True

        # igore 2nd law for biomass transport
        biomass_present = [met.biomass for met in self.metabolites]
        if all(biomass_present) if bool(biomass_present) else False:
            self._ignore_snd = True
        
        for met, stoich in self.metabolites.items(): #initialisation will be slightly longer but it avoids to compute it each time net_elements is called
            cpd=drg_tools.get_compound(met)
            #if self._cobra_formula==False:
             #   self._cobra_formula=True if (cpd is None) or (cpd.can_be_transformed() is False) or (cpd.atom_bag is None) or (cpd.atom_bag == {}) else False

    def net_elements(self, balance_mg = False, round_dp = False,rxn_already_balanced = True):
        return drg_tools.net_elements(self, balance_mg = balance_mg, round_dp = round_dp,rxn_already_balanced = rxn_already_balanced)
    
    def transported_c_h(self,  round_dp = False):
        return drg_tools.transported_c_h(self,  round_dp = round_dp)
    
    def _cobra_formula(self):
        """ is used in several methods of drG_tools : if one metabolite needs to use cobra formula instead of equilibrator, 
        ###calc_average_charge_proton will be based on all metabolites' cobra formula"""            
        return self._cobra_formula
    

    def check_atom_balance(self, round_dp = False, electrons = False, pMg=Q_(14,'')):
        '''Check the mass balance of a reaction 

        Uses the eQuilibrator compound atom bag to check the mass balance. If this is not available the COBRA model formula is used. 
        Note magnesium balance is ignored by default.
        
        Parameters
        ----------
        rxn : cobra.Reaction
            Reaction to check
        round_dp : int
            Number of decimal places to round to
        electrons : bool
            Include electrons in the balance
        pMg : float
            pMg of the reaction
        
        Returns
        -------
        net_elements : dict
             a dict of {element: amount} for unbalanced elements.
            This should be empty for balanced reactions.
            
        '''
        rxn=self
        atom_bags = {}

        for met, stoich in rxn.metabolites.items():
            if met.compound is not None:
                if met.compound.atom_bag is not None:
                    atom_bag = {element:value*stoich for element, value in met.compound.atom_bag.items()}
                    atom_bag['H'] = met.average_charge_protons(round_dp = round_dp, pMg=pMg)[1]*stoich                
                    atom_bags[met.id] = atom_bag 

                else:
                    atom_bags[met.id] = {element:value*stoich for element, value in met.elements.items()}

            else:
                atom_bags[met.id] = {element:value*stoich for element, value in met.elements.items()}

        rxn_atom_bags = pd.DataFrame.from_dict(atom_bags, orient="columns").fillna(0)

        net_elements = rxn_atom_bags.sum(axis=1)

        #remove 0 values
        net_elements = net_elements[net_elements != 0].to_dict()
        #remove values with abs < 1e-6
        net_elements = {key: value for key, value in net_elements.items() if abs(value) > 1e-6}
        if not electrons:
            net_elements.pop('e-', None)

        return net_elements

    @property
    def ignore_snd(self):
        """Flag to ignore the second law constraint for a reaction"""
        return self._ignore_snd

    @ignore_snd.setter
    def ignore_snd(self, value):
        self._ignore_snd = value


    @property
    def drG0(self):
        """Standard Gibbs energy of a reaction before transformation."""
        return self._drG0

    @drG0.setter
    def drG0(self, value):
        self._drG0 = value

    @property
    def drG0prime(self):
        """Transformed Gibbs energy of a reaction."""
        return self._drG0prime

    @drG0prime.setter
    def drG0prime(self, value):
        self._drG0prime = value

    @property
    def drGmprime(self):
        """Transformed Gibbs energy of a reaction."""
        return self._drGmprime

    @drGmprime.setter
    def drGmprime(self, value):
        self._drGmprime = value

    @property
    def drGtransport(self):
        """Transport componenent of gibbs energy of a transporter reaction."""
        return self._drGtransport

    @drGtransport.setter
    def drGtransport(self, value):
        self._drGtransport = value

    @property
    def drG_h_transport(self):
        """Proton transport component of the transport component of the gibbs energy of a transporter reaction"""
        return self._drG_h_transport

    @drG_h_transport.setter
    def drG_h_transport(self, value):
        self._drG_h_transport = value

    @property
    def drG_c_transport(self):
        """Charge transport component of the transport component of the gibbs energy of a transporter reaction"""
        return self._drG_c_transport

    @drG_c_transport.setter
    def drG_c_transport(self, value):
        self._drG_c_transport = value

    @property
    def drGtransform(self):
        """Gibbs energy to be added to the standard gibbs energy of a reaction to get the transformed gibbs free energy"""

        return self._drGtransform

    @drGtransform.setter
    def drGtransform(self, value):
        self._drGtransform = value

    @property
    def drG(self):
        """Gibbs free energy of the reaction. drG0' + RTlnC"""
        return self._drG

    @drG.setter
    def drG(self, value):
        self._drG = value
    
    @property
    def drG_SE(self):
        """standard error on the gibbs free energy estimate. Can be set explicitly or calcualted from the covariance"""
        return self._drG_SE

    @drG_SE.setter
    def drG_SE(self, value):
        self._drG_SE = value

    @property
    def transported_h(self):
        """Free or additional protons transported in a transport reaction"""
        return self._transported_h

    @transported_h.setter
    def transported_h(self, value, round_dp = False):
        #check correct compartments have been defined 
        if value.keys() == set(self.compartments):

            #only balance charge and mg if they are already in the reaction 
            balance_charge = False
            if len([met for met in self.metabolites if (met in self.model.charge_dict.values())]) != 0:
                balance_charge = True

            balance_mg = False
            if len([met for met in self.metabolites if (met in self.model.mg_dict.values())]) != 0:
                balance_mg = True

            i = 0 
            for stoich in value.values():
                i+= stoich

            if i != 0:
                corrected_dict = {comp:((i*-2)+1)*list(value.values())[0] for i, comp in enumerate(value.keys())}
                warnings.warn(f"Transported free protons do not balance. Corrected to {corrected_dict}", stacklevel=2)
                self._transported_h = corrected_dict

                drg_tools.reaction_balance(self, balance_charge=balance_charge, balance_mg=balance_mg, round_dp = round_dp)

            else:
                self._transported_h = value
                drg_tools.reaction_balance(self, balance_charge=balance_charge, balance_mg=balance_mg, round_dp = round_dp)

        else:
            warnings.warn(f"Incorrect compartments defined. Reaction compartments are {self.compartments}", stacklevel=2)

    @property
    def transported_charge(self):
        """Free or additional charge transported in a transport reaction. 
        Represents a generic positive ion"""

        #if the reaction is a generic charge transporter transported_charge and the reaction sotichimetry must match
        if  all([met in self.model.charge_dict.values()
                    for met in self.metabolites]):
            
            self._transported_charge = {met.compartment:stoich for met, stoich in self.metabolites.items()}

        return self._transported_charge

    @transported_charge.setter
    def transported_charge(self, value):
        #check correct compartments have been defined 
        if value.keys() == set(self.compartments):

            #only balance charge and mg if they are already in the reaction 
            balance_charge = False
            if len([met for met in self.metabolites if (met in self.model.charge_dict.values())]) != 0:
                balance_charge = True

            balance_mg = False
            if len([met for met in self.metabolites if (met in self.model.mg_dict.values())]) != 0:
                balance_mg = True

            i = 0 
            for stoich in value.values():
                i+= stoich

            if i != 0:
                if not self.boundary:  # ignore this warning for boundary reactions
                    corrected_dict = {comp:((i*-2)+1)*list(value.values())[0] for i, comp in enumerate(value.keys())}
                    warnings.warn(f"Transported free charge does not balance. Corrected to {corrected_dict}", stacklevel=2)
                    self._transported_charge = corrected_dict

                    drg_tools.reaction_balance(self, balance_charge=balance_charge, balance_mg=balance_mg)

            else:
                self._transported_charge = value
                drg_tools.reaction_balance(self, balance_charge=balance_charge, balance_mg=balance_mg)

        else:
            warnings.warn(f"Incorrect compartments defined. Reaction compartments are {self.compartments}", stacklevel=2)

    @property
    def transported_mets(self):
        """explicitly define metabolite transported in a reaction. This is useful for
        transporters that include chemical transformation during transport. Leave empty for 
        automatic calculation of transported metabolites"""
        return self._transported_mets

    @transported_mets.setter
    def transported_mets(self, value):
        self._transported_mets = value

    @property
    def balanced(self):
        """Defines if a reaction has already been charge and proton balanced"""
        return self._balanced

    @balanced.setter
    def balanced(self, value):
        self._balanced = value

    def split_reaction(self):
        """Split a multiple compartment reaction into subunits of 2 compartments or less"""
        rxn_info=pd.DataFrame(columns=['name','compartment','stoichiometric coefficient'])
        for met,stoich in self.metabolites.items():
            #we want the name of the metabolite without the compartment and to store relevant info in a df for later
            name='_'.join(met.id.split('_')[:-1])
            rxn_info.loc[met.id,'name']=name
            rxn_info.loc[met.id,'compartment']=met.compartment
            rxn_info.loc[met.id,'stoichiometric coefficient']=stoich
        #check for duplicates in name column == metabolites that are present in 2 compartments
        #and build sub reaction with conserved metabolites over transport of 2 compartments
        list_sub_rxn=[]
        for i,row in rxn_info[rxn_info.duplicated(subset=['name'])]['name'].items():
            #we get the name of the metabolite that is present twice 
            sub_rxn=ThermoReaction(Reaction('sub_rxn_transport_'+str(row)))
            sub_rxn.add_metabolites({self.model.metabolites.get_by_id(rxn_info[rxn_info['name']==row].index[0]):rxn_info[rxn_info['name']==row]['stoichiometric coefficient'].values[0],self.model.metabolites.get_by_id(rxn_info[rxn_info['name']==row].index[1]):rxn_info[rxn_info['name']==row]['stoichiometric coefficient'].values[1]})
            list_sub_rxn.append(sub_rxn)
            #remove the metabolites from the original reaction
            self.subtract_metabolites({self.model.metabolites.get_by_id(rxn_info[rxn_info['name']==row].index[0]):rxn_info[rxn_info['name']==row]['stoichiometric coefficient'].values[0],self.model.metabolites.get_by_id(rxn_info[rxn_info['name']==row].index[1]):rxn_info[rxn_info['name']==row]['stoichiometric coefficient'].values[1]})
        #We now have : sub reactions with conserved metabolites and the initial reaction without the conserved metabolites
        #we balance the sub reactions and the modified initial reaction
        #then we add the balanced sub reactions to the modified initial reaction
        return list_sub_rxn