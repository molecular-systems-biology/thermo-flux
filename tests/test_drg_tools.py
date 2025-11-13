import pytest
import numpy as np
from thermo_flux.tools import drg_tools
from equilibrator_api import ComponentContribution, Q_

cc = ComponentContribution()

eps = 1e-5

class TestDrgTools:
    
    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_get_compound(self, tmodel):
        cpd = drg_tools.get_compound(tmodel.metabolites[1])
        assert cpd == cc.get_compound('Kegg:C00002')

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_dfG0(self, tmodel):
        dfG0 = tmodel.metabolites[1].dfG0
        print(dfG0)
        test_dfG0 = Q_(-2811.578331958078, "kJ/mol")
        assert abs(dfG0.m - test_dfG0.m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_calc_dfG_transform(self, tmodel):
        dfG_transform = drg_tools.calc_dfG_transform(tmodel.metabolites[1])
        test_dfG_transform = Q_(521.9624321661707, "kJ/mol")
        assert abs(dfG_transform.m - test_dfG_transform.m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_dfG0prime(self, tmodel):
        dfG0prime = tmodel.metabolites[1].dfG0prime
        test_dfG0prime = Q_(-2289.615899791907, "kJ/mol")
        assert abs(dfG0prime.m - test_dfG0prime.m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_dfG0_cov_sqrt(self, tmodel):
        """test that metabolites in different compartments are correlated"""
        id_list = [met.id[:-2] for met in tmodel.metabolites if met.id[:-2] not in ['charge', 'Mg']]
        multiples = set(a for a in set(id_list) if id_list.count(a) > 1)

        multi_mets = []
        for met in multiples:
            matches = [i for i, e in enumerate(id_list) if e == met]
            multi_mets.append(matches)

        correlated = []
        for pair in multi_mets:
            correlated.append(np.isclose(tmodel.dfG0_cov_sqrt[pair], tmodel.dfG0_cov_sqrt[pair[0]]).all())
        
        assert all(correlated)

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_drG0(self, tmodel):
        drG0 = tmodel.reactions[0].drG0
        test_drG0 = Q_(52.07649473101105, "kJ/mol")
        assert abs(drG0.m - test_drG0.m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_drG0prime(self, tmodel):
        drG0prime = tmodel.reactions[0].drG0prime
        test_drG0prime = Q_(-35.72304430928557, "kJ/mol")
        assert abs(drG0prime.m - test_drG0prime.m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_transported_c_h(self, tmodel):
        h_c = drg_tools.transported_c_h(tmodel.reactions.get_by_id("PYRt-1"))
        assert h_c == (3.0, -1.0, -3.0, 1.0, 'c', 'e', True)

    #test inner compartment detection
    def test_inner_compartment_detection(self, tmodel):
        rxn = tmodel.reactions.get_by_id("PYRt-1")
        inner_comp, outer_comp = drg_tools._transport_direction(rxn)
        assert inner_comp == 'c'
        assert outer_comp == 'e'
        assert rxn.inner_compartment == 'c'

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_drG_transport(self, tmodel):
        drG_transport, drg_h, drg_c = drg_tools.calc_drGtransport(tmodel.reactions.get_by_id("PYRt-1"))
        assert abs(drG_transport.m - Q_(-29.014609533125487, "kJ/mol").m) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_net_elements(self, tmodel):
        net_elements, transported_free_h = drg_tools.net_elements(tmodel.reactions.ASNS1)
        
        expected = {tmodel.metabolites.h_c: -2,
                   tmodel.metabolites.charge_c: -1.6135891238290512,
                   tmodel.metabolites.Mg_c: 0.1932054380854744}
        
        for met, val in expected.items():
            assert abs(net_elements[met] - val) < eps

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_net_elements_transport(self, tmodel):
        net_elements, transported_free_h = drg_tools.net_elements(tmodel.reactions.ASNS1)

        assert net_elements == {tmodel.metabolites.h_c: -2,
                               tmodel.metabolites.charge_c: -1.6135891238290512,
                               tmodel.metabolites.Mg_c: 0.1932054380854744}

    @pytest.mark.usefixtures("tmodel")  # Use the default model
    def test_reaction_balance(self, tmodel):
        drg_tools.reaction_balance(tmodel.reactions.ASNS1, balance_mg=False)

        assert tmodel.reactions.ASNS1.metabolites == {
            tmodel.metabolites.get_by_id("asp-L_c"): -1.0,
            tmodel.metabolites.atp_c: -1.0,
            tmodel.metabolites.get_by_id("gln-L_c"): -1.0,
            tmodel.metabolites.h2o_c: -1.0,
            tmodel.metabolites.amp_c: 1.0,
            tmodel.metabolites.get_by_id("asn-L_c"): 1.0,
            tmodel.metabolites.get_by_id("glu-L_c"): 1.0,
            tmodel.metabolites.ppi_c: 1.0,
            tmodel.metabolites.h_c: 1.687350491704554
        }

    @pytest.mark.usefixtures("tmodel_b")  # Use the alternative model
    def test_ATPs(self, tmodel_b):
        rxn = tmodel_b.reactions.get_by_id('ATPS3m')

        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        test_balanced_mets = {
            tmodel_b.metabolites.adp_m: -3.0,
            tmodel_b.metabolites.pi_m: -3.0,
            tmodel_b.metabolites.atp_m: 3.0,
            tmodel_b.metabolites.h2o_m: 3.0,
            tmodel_b.metabolites.h_c: -10.0,
            tmodel_b.metabolites.h_m: 11.555101,
            tmodel_b.metabolites.charge_c: -10.0,
            tmodel_b.metabolites.charge_m: 10.0
        }

        assert rxn.metabolites == test_balanced_mets

        rxn.transported_h = {'c': -20.0, 'm': 20.0}

        test_balanced_mets = {
            tmodel_b.metabolites.adp_m: -3.0,
            tmodel_b.metabolites.pi_m: -3.0,
            tmodel_b.metabolites.atp_m: 3.0,
            tmodel_b.metabolites.h2o_m: 3.0,
            tmodel_b.metabolites.h_c: -20.0,
            tmodel_b.metabolites.h_m: 21.555101,
            tmodel_b.metabolites.charge_c: -20.0,
            tmodel_b.metabolites.charge_m: 20.0
        }

        assert rxn.metabolites == test_balanced_mets

        # test rebalancing doesn't cause issues
        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        assert rxn.metabolites == test_balanced_mets

    @pytest.mark.usefixtures("tmodel_b")  # Use the alternative model
    def test_PIt_1(self, tmodel_b):
        '''PI crosses the membrane, but the major metabolite in the inner compartment (c) is HPO4^2-'''
        rxn = tmodel_b.reactions.get_by_id('PIt-1')
        drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        test_balanced_mets = {
            tmodel_b.metabolites.pi_e: -1.0,
            tmodel_b.metabolites.charge_c: -1.0,
            tmodel_b.metabolites.pi_c: 1.0,
            tmodel_b.metabolites.charge_e: 1.0,
            tmodel_b.metabolites.h_c: 0.845221
        }

        assert rxn.metabolites == test_balanced_mets
