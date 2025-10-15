import unittest
from cobra.io import read_sbml_model
from thermo_flux.core.model import ThermoModel
from thermo_flux.io import load_excel as ex
from equilibrator_api import Q_


class TestModel(unittest.TestCase):

    maxDiff = None

    @classmethod
    def setUpClass(cls):

        cls.model = read_sbml_model('yeast_test.xml')
        cls.tmodel = ThermoModel(cls.model,
                                 pH={'c': Q_(7),
                                     'm': Q_(7.4),
                                     'e': Q_(5)},
                                 I={'c': Q_(0.25, 'M'),
                                    'm': Q_(0.25, 'M'),
                                    'e': Q_(0.25, 'M')},
                                 T=Q_(303.15, 'K'),
                                 pMg={'c': Q_(3), 'm': Q_(3), 'e': Q_(3)},
                                 phi={'ec': Q_(-0.06, 'V'),
                                      'cm': Q_(-0.16, 'V')},
                                 update_thermo_info=False)

    def test_get_rxn(self):
        rxn = (self.tmodel.reactions[0])
        self.assertEqual(rxn.id, 'ASNS1')

        return 

    def test_model_thermo_params(self):
        self.assertTrue(hasattr(self.tmodel, "pH"))
        self.assertTrue(hasattr(self.tmodel, "I"))
        self.assertTrue(hasattr(self.tmodel, "pMg"))
        self.assertTrue(hasattr(self.tmodel, "T"))
        self.assertTrue(hasattr(self.tmodel, "phi"))
        self.assertTrue(hasattr(self.tmodel, "gdiss_lim"))
        
        return

    def test_phi_dict(self):
        phi_dict = self.tmodel.phi_dict
        test_phi_dict = {'c': {'c': Q_(0,'V'), 'm': Q_(-0.16,'V'), 'e': Q_(0.06,'V')},
                         'm': {'c': Q_(0.16,'V'), 'm': Q_(0,'V'), 'e': Q_(0,'V')},
                         'e': {'c': Q_(-0.06,'V'), 'm': Q_(0,'V'), 'e': Q_(0,'V')}}

        self.assertEqual(phi_dict, test_phi_dict)
        return

    def test_rxn_thermo_params(self):
        for rxn in self.tmodel.reactions:
            self.assertTrue(hasattr(rxn, "drG0"))
            self.assertTrue(hasattr(rxn, "drG0prime"))
            self.assertTrue(hasattr(rxn, "drGtransport"))
            self.assertTrue(hasattr(rxn, "drGtransform"))
            self.assertTrue(hasattr(rxn, "drG_h_transport"))
            self.assertTrue(hasattr(rxn, "drG_c_transport"))
            self.assertTrue(hasattr(rxn, "drG"))
            self.assertTrue(hasattr(rxn, "drG_SE"))          

        return 

    def test_met_thermo_params(self):
        for met in self.tmodel.metabolites:
            self.assertTrue(hasattr(met, "upper_bound"))
            self.assertTrue(hasattr(met, "lower_bound"))
            self.assertTrue(hasattr(met, "concentration"))
            self.assertTrue(hasattr(met, "accession"))
            self.assertTrue(hasattr(met, "dfG0"))
            self.assertTrue(hasattr(met, "dfG0prime"))       

        return 

    def test_proton_dict(self):
        proton_dict = self.tmodel.proton_dict
        self.assertEqual(len(proton_dict), len(self.tmodel.compartments))
        return

    def test_charge_dict(self):
        charge_dict = self.tmodel.charge_dict
        self.assertEqual(len(charge_dict), len(self.tmodel.compartments))
        return

    def test_mg_dict(self):
        mg_dict = self.tmodel.mg_dict
        self.assertEqual(len(mg_dict), len(self.tmodel.compartments))
        return


class TestExcelModel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = ex.create_cobra_model('ecoli', 
                                          model_excel="..\\examples\\ecoli\\model.xlsx", 
                                          keggids_csv="..\\examples\\ecoli\\ecoli_kegg_id.csv")

        cls.tmodel = ThermoModel(cls.model)

    def test_pH(self):
        pH = self.tmodel.pH
        self.assertEqual(pH, {'c': Q_(7.6), 'e': Q_(7)})

    def test_T(self):
        T = self.tmodel.T
        self.assertEqual(T, Q_(310.15, 'K'))

    def test_phi_dict(self):
        phi_dict = self.tmodel.phi_dict
        test_phi_dict = {'c': {'c': Q_(0, 'V'), 'e': Q_(-0.15, 'V')},
                         'e': {'c': Q_(0.15, 'V'), 'e': Q_(0, 'V')}}
        self.assertEqual(phi_dict, test_phi_dict)

if __name__ == '__main__':
    unittest.main()
