"""A metabolite with thermodynamic properties"""

from cobra.core import Metabolite

# from thermo_flux.utils.annotation import get_suitable_ids

from thermo_flux.tools import drg_tools 

from equilibrator_api import Q_
from equilibrator_cache.models.compound import Compound as eQ_compound


from warnings import warn



proton_annotations = {'Kegg': 'C00080',
     'bigg.metabolite': 'h',
     'chebi': 'CHEBI:15378',
     'hmdb': 'HMDB59597',
     'kegg': 'C00080',
     'metacyc.compound': 'PROTON',
     'reactome': '1132304',
     'sabiork.compound': '39',
     'seed': 'cpd00067',
     'metanetx.chemical': 'MNXM145872',
     'synonyms': 'H+',
     'biocyc': 'META:PROTON',
     'inchi_key': 'GPRLSGONYQIRFK-UHFFFAOYSA-N',
     'metanetx.chemical': 'MNXM1',
     'seed.compound': 'cpd00067'}

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


class ThermoMetabolite(Metabolite):
    """Thermodynamic metabolite object.

    takes cobra metabolite and adds thermodynamic properties

    """

    # add concentration variable if using optlang 

    def __init__(self, metabolite, model=None,
                 upper_bound: Q_ = Q_(10, 'mM'),
                 lower_bound: Q_ = Q_(0.1, 'uM'), #ED leave default lower bound as 0.1 uM (this is already quite a low concnetration! check what happens for auto analysis of all models...)
                 concentration: Q_ = Q_(1, 'M'), 
                 accession: str = None,
                 dfG0: Q_ = None,
                 dfG0prime: Q_ = None,
                 compound=None,
                 dfG_SE: Q_ = None,
                 redox: bool = False,
                 biomass: bool = False,
                 unknown: bool = False,
                 ignore_conc: bool = False):

        Metabolite.__init__(self)
        self._concentration = concentration
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._model = model
        self._accession = accession
        self._dfG0 = dfG0
        self._dfG0prime = dfG0prime
        self._compound = compound
        self._dfG_SE = dfG_SE
        self._redox = redox
        self._biomass = biomass
        self._unknown = unknown 
        self._ignore_conc = ignore_conc
        for attr, value in metabolite.__dict__.items():
            if attr not in ['_model', '_reaction']:
                self.__dict__[attr] = value

        if self.annotation != {}:
            if any([proton_annotations[identifier] in key if isinstance(key, list) else key == proton_annotations[identifier] 
                      for identifier, key in self.annotation.items()
                      if identifier in proton_annotations]):
                self._lower_bound = Q_(1, 'M')
                self._upper_bound = Q_(1, 'M')
                self._ignore_conc = True

            if any([water_annotations[identifier] in key if isinstance(key, list) else key == water_annotations[identifier] 
                      for identifier, key in self.annotation.items()
                      if identifier in water_annotations]):
                
                self._lower_bound = Q_(1, 'M')
                self._upper_bound = Q_(1, 'M')
                self._ignore_conc = True  

        if 'biomass' in self.id.lower() or self.biomass:
            self._biomass = True
            self._ignore_conc = True

    def __eq__(self, other):
        '''sometimes metabolite are copied and then have a different memory address. override eq to prevent issues'''
        return isinstance(other, ThermoMetabolite) and self.id == other.id

    def __hash__(self):
        return hash(self.id)

    @property
    def lower_bound(self):
        return self._lower_bound 

    @lower_bound.setter
    def lower_bound(self, value):
        self._lower_bound = value
        self._check_bounds(value, self._upper_bound)
        
    @property
    def upper_bound(self):
        return self._upper_bound
        
    @upper_bound.setter
    def upper_bound(self, value):
        self._upper_bound = value

    @property
    def concentration(self):
        return self._concentration
        
    @concentration.setter
    def concentration(self, value):
        self._concentration = value

    @staticmethod
    def _check_bounds(lb, ub):
        if ub < lb:
            raise ValueError(
                "The lower bound must be less than or equal to the upper bound ({} <= {}).".format(lb, ub))

    @property
    def accession(self):
        if self._accession is None:
            if self.compound is not None:
                self._accession  = self.compound.id
                #self._accession = self.id[:-2] #ToDo need faster way of updating accessions 
            else:
                self._accession = self.id[:-2] #ToDo need more robust way of defining unknown metabolites in a compartment agnostic way
        return self._accession

    @accession.setter
    def accession(self, value):
        self._accession = value

    @property
    def dfG0(self):
        return self._dfG0

    @dfG0.setter
    def dfG0(self, value):
        self._dfG0 = value

    @property
    def dfG0prime(self):
        return self._dfG0prime

    @dfG0prime.setter
    def dfG0prime(self, value):
        self._dfG0prime = value

    @property
    def dfGprime(self):
        warn('dfGprime is deprecated use dfG0prime instead', DeprecationWarning, stacklevel=2)
        return self._dfG0prime

    @dfGprime.setter
    def dfGprime(self, value):
        warn('dfGprime is deprecated use dfG0prime instead', DeprecationWarning, stacklevel=2)
        self._dfG0prime = value

    @property
    def compound(self, update_annotations = False):
        if self._compound is None:
            cpd, annotation, formula, inchi, searched = drg_tools.get_suitable_ids(self, update_annotations = update_annotations)
            self.annotation.update(annotation)
            self._compound = cpd
        return self._compound

    @compound.setter #only accept eQ_compound objects or None
    def compound(self, value):
        if not isinstance(value, eQ_compound) and value is not None:
            raise ValueError('compound must be an instance of equilibrator_cache.models.compound.Compound')
        self._compound = value


    @property
    def dfG_SE(self):
        return self._dfG_SE

    @dfG_SE.setter
    def dfG_SE(self, value):
        #dfG_SE must be the same for multiple occurances of a metabolite in the model 
        #if it is updated in one place it must automatically update in the other
        if self.model is not None:
            for met in self.model.metabolites:
                if met.accession == self.accession: 
                    met._dfG_SE = value
        self._dfG_SE = value


    @property
    def redox(self):
        return self._redox

    @redox.setter
    def redox(self, value):
        self._redox = value

    @property
    def biomass(self):
        return self._biomass

    @biomass.setter
    def biomass(self, value):
        self._biomass = value

    @property
    def unknown(self):
        return self._unknown

    @unknown.setter
    def unknown(self, value):
        self._unknown = value
        if value is True:
            self._compound = None

    @property
    def ignore_conc(self):
        return self._ignore_conc

    @ignore_conc.setter
    def ignore_conc(self, value):
        self._ignore_conc = value

    def average_charge_protons(self,  pH = None, pMg=None, ionic_strength=None, temperature=None, accuracy = 0.1, round_dp=False,cobra_formula=False):
        if pH is None:
            pH = self.model.pH[self.compartment]

        if pMg is None:
            pMg = self.model.pMg[self.compartment]

        if ionic_strength is None:
            ionic_strength = self.model.I[self.compartment]

        if temperature is None:
            temperature = self.model.T

        return drg_tools.calc_average_charge_protons(self, pH, pMg, ionic_strength, temperature, accuracy, round_dp,cobra_formula=cobra_formula)

    def major_microspecies(self, pH=None, pMg=None, ionic_strength=None, temperature=None):
        if pH is None:
            pH = self.model.pH[self.compartment]

        if pMg is None:
            pMg = self.model.pMg[self.compartment]

        if ionic_strength is None:
            ionic_strength = self.model.I[self.compartment]

        if temperature is None:
            temperature = self.model.T

        return drg_tools.major_microspecies(self, pH, pMg, ionic_strength, temperature)


    def check_consistency(self, ignore_H = True):
        '''Check the consistency of a metabolite between the original model definiton and the metabolite identified in equilibrator database for thermodynamic analysis
        '''
        eq_atom_bag = {}
        if self.compound is not None:
            if self.compound.atom_bag is not None:
                eq_atom_bag = self.compound.atom_bag
                eq_atom_bag.pop('e-', None)  
                if ignore_H:
                    eq_atom_bag.pop('H', None)
            else: #if the equilibrator compound does not have an atom bag then return True regardelss of cobra model formula 
                return True #this ensures some metabolites with ambiguous structures that are in the equilibrator databse are still correctly assigned 
        
        cobra_elements = {}
        if self.elements is not None:
            cobra_elements = self.elements
        if ignore_H:
            cobra_elements.pop('H', None)

        if cobra_elements == {}:
            return True #if no cobra formula is set then defult to True and equilibrator formula will be used from then on

        if cobra_elements != eq_atom_bag:
            return False
        else:
            return True
