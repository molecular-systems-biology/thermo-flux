import unittest
import numpy as np
from cobra.io import read_sbml_model
import sys
sys.path.append('/home/nathan/public_thermo_flux/')
from thermo_flux.core.model import ThermoModel
from thermo_flux.tools import drg_tools
from equilibrator_api import ComponentContribution, Q_

cc = ComponentContribution()

eps=1e-5

class TestModel(unittest.TestCase):
    maxDiff = None


    @classmethod
    def setUpClass(cls):

        cls.model = read_sbml_model('yeast_merged.xml')
        cls.tmodel = ThermoModel(cls.model, 
                                 pH={'c': Q_(7), 'm': Q_(7.4), 'e': Q_(5)},
                                 I={'c': Q_(0.25, 'M'), 'm': Q_(0.25, 'M'), 
                                    'e': Q_(0.25, 'M')}, 
                                 T=Q_(303.15, 'K'), 
                                 pMg={'c': Q_(3), 'm': Q_(3), 'e': Q_(3)},
                                 phi={'ec': Q_(-0.06, 'V'),
                                      'cm': Q_(-0.16, 'V')}, 
                                 update_thermo_info=True)

        cls.model = read_sbml_model('yeast_merged.xml')
        cls.tmodel_b = ThermoModel(cls.model,pH = {'c':Q_(7), 'm':Q_(1),'e':Q_(5)},
                                            I = {'c':Q_(0.25,'M'), 'm':Q_(0.25,'M'),'e':Q_(0.25,'M')},
                                            T = Q_(303.15,'K') ,
                                            pMg = {'c':Q_(3), 'm':Q_(3),'e':Q_(3)},
                                            phi = {'ec':Q_(-0.06,'V'), 'cm':Q_(-0.16,'V')},
                                            update_thermo_info = False)


    def test_get_compound(self):
        cpd = drg_tools.get_compound(self.tmodel.metabolites[1])
        self.assertEqual(cpd, cc.get_compound('Kegg:C00002'))


    def test_dfG0(self):
        dfG0 = self.tmodel.metabolites[1].dfG0
        test_dfG0 = Q_(-2811.578331958078,"kJ/mol")

        self.assertAlmostEqual(dfG0, test_dfG0,delta=eps)

    def test_calc_dfG_transform(self):  # noqa: W191

        dfG_transform = drg_tools.calc_dfG_transform(self.tmodel.metabolites[1])
        test_dfG_transform = Q_(521.9624321661707,"kJ/mol")

        self.assertAlmostEqual(dfG_transform, test_dfG_transform,delta=eps)

    def test_dfG0prime(self):
        dfG0prime = self.tmodel.metabolites[1].dfG0prime
        test_dfG0prime = Q_(-2289.615899791907,"kJ/mol")

        self.assertAlmostEqual(dfG0prime, test_dfG0prime,delta=eps)

    def test_dfG0_cov_sqrt(self):
        """test that metabolites in different comparmtents are correlated"""

        id_list = [met.id[:-2] for met in self.tmodel.metabolites if met.id[:-2] not in ['charge', 'Mg']] #use met id to find matches
        multiples = set(a for a in set(id_list) if id_list.count(a) > 1) #all metabolites that appear more than once

        multi_mets = [] #list of lists of met indexes that occur in multiple compartments 
        for met in multiples:
            matches = [i for i, e in enumerate(id_list) if e == met]
            multi_mets.append(matches)

        correlated = [] 
        for pair in multi_mets:
            correlated.append(np.isclose(self.tmodel.dfG0_cov_sqrt[pair], self.tmodel.dfG0_cov_sqrt[pair[0]]).all())
        
        self.assertTrue(all(correlated))

    def test_drG0(self):
        drG0 = self.tmodel.reactions[0].drG0
        test_drG0 = Q_(52.07649473101105,"kJ/mol")
        self.assertAlmostEqual(drG0.m, test_drG0.m,delta=eps)

    def test_drG0prime(self):
        drG0prime = self.tmodel.reactions[0].drG0prime
        test_drG0prime = Q_(-35.72304430928557,"kJ/mol")
        self.assertAlmostEqual(drG0prime.m, test_drG0prime.m,delta=eps)

    def test_transported_c_h(self):
        h_c = drg_tools.transported_c_h(self.tmodel.reactions.get_by_id("PYRt-1"))
        self.assertEqual(h_c, (3.0, -1.0, -3.0, 1.0, 'c', 'e', True))

    def test_drG_transport(self):
        drG_transport, drg_h, drg_c = drg_tools.calc_drGtransport(self.tmodel.reactions.get_by_id("PYRt-1"))
        self.assertAlmostEqual(drG_transport.m, Q_(-29.014609533125487, "kJ/mol").m,delta=eps)

    def test_net_elements(self):
        net_elements, transported_free_h = drg_tools.net_elements(self.tmodel.reactions.ASNS1)
        
        self.assertAlmostEqual(net_elements, {self.tmodel.metabolites.h_c:-2,
                                        self.tmodel.metabolites.charge_c:-1.6135891238290512,
                                        self.tmodel.metabolites.Mg_c:0.1932054380854744})

    def test_net_elements_transport(self):
        net_elements, transported_free_h = drg_tools.net_elements(self.tmodel.reactions.ASNS1)

        self.assertEqual(net_elements, {self.tmodel.metabolites.h_c:-2,
                                        self.tmodel.metabolites.charge_c:-1.6135891238290512,
                                        self.tmodel.metabolites.Mg_c:0.1932054380854744})

    def test_reaction_balance(self):
        drg_tools.reaction_balance(self.tmodel.reactions.ASNS1, balance_mg=False)

        self.assertEqual(self.tmodel.reactions.ASNS1.metabolites,
                         {self.tmodel.metabolites.get_by_id("asp-L_c"): -1.0,
                          self.tmodel.metabolites.atp_c: -1.0,
                          self.tmodel.metabolites.get_by_id("gln-L_c"): -1.0,
                          self.tmodel.metabolites.h2o_c: -1.0,
                          self.tmodel.metabolites.amp_c: 1.0,
                          self.tmodel.metabolites.get_by_id("asn-L_c"): 1.0,
                          self.tmodel.metabolites.get_by_id("glu-L_c"): 1.0,
                          self.tmodel.metabolites.ppi_c: 1.0,
                          self.tmodel.metabolites.h_c: 1.687350491704554})
                        
    # test reaction coupled transporter
    def test_ATPs(self):
        rxn = self.tmodel_b.reactions.get_by_id('ATPS3m')

        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        test_balanced_mets = {self.tmodel_b.metabolites.adp_m: -3.0,
                              self.tmodel_b.metabolites.pi_m: -3.0,
                              self.tmodel_b.metabolites.atp_m: 3.0,
                              self.tmodel_b.metabolites.h2o_m: 3.0,
                              self.tmodel_b.metabolites.h_c: -10.0,
                              self.tmodel_b.metabolites.h_m: 11.555101,
                              self.tmodel_b.metabolites.charge_c: -10.0,
                              self.tmodel_b.metabolites.charge_m: 10.0}

        self.assertEqual(rxn.metabolites, test_balanced_mets)

        rxn.transported_h = {'c': -20.0, 'm': 20.0}

        test_balanced_mets = {self.tmodel_b.metabolites.adp_m: -3.0,
                              self.tmodel_b.metabolites.pi_m: -3.0,
                              self.tmodel_b.metabolites.atp_m: 3.0,
                              self.tmodel_b.metabolites.h2o_m: 3.0,
                              self.tmodel_b.metabolites.h_c: -20.0,
                              self.tmodel_b.metabolites.h_m: 21.555101,
                              self.tmodel_b.metabolites.charge_c: -20.0,
                              self.tmodel_b.metabolites.charge_m: 20.0}

        self.assertEqual(rxn.metabolites, test_balanced_mets)
 

        # test rebalancing doesn't cause issues 
        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        self.assertEqual(rxn.metabolites, test_balanced_mets)

    # test transporter with additional protons 
    def test_PIt_1(self):
        '''PI crosses the membrane, but the major metabolite in the inner compartment (c) is HPO4^2- and is transported according to our assumptions.
          Thus this transport from e to c leaves 1 protons in e and needs 0.2 additional protons from c, balancing with the proton symport h_e -> h_c  '''
        rxn = self.tmodel_b.reactions.get_by_id('PIt-1')
        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        test_balanced_mets = {self.tmodel_b.metabolites.pi_e: -1.0,
                              self.tmodel_b.metabolites.charge_c: -1.0,
                              self.tmodel_b.metabolites.pi_c: 1.0,
                              self.tmodel_b.metabolites.charge_e: 1.0,
                            #   self.tmodel_b.metabolites.h_e: .0,
                              self.tmodel_b.metabolites.h_c: 0.845221} ##actually

        self.assertEqual(rxn.metabolites, test_balanced_mets)
if __name__ == '__main__':
    unittest.main()
