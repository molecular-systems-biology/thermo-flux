import pytest
from cobra.io import read_sbml_model
from thermo_flux.core.model import ThermoModel
from thermo_flux.io import load_excel as ex
from equilibrator_api import Q_


class TestModel:

    @pytest.mark.usefixtures("tmodel")
    def test_get_rxn(self, tmodel):
        rxn = tmodel.reactions[0]
        assert rxn.id == 'ASNS1'

    @pytest.mark.usefixtures("tmodel")
    def test_model_thermo_params(self, tmodel):
        assert hasattr(tmodel, "pH")
        assert hasattr(tmodel, "I")
        assert hasattr(tmodel, "pMg")
        assert hasattr(tmodel, "T")
        assert hasattr(tmodel, "phi")
        assert hasattr(tmodel, "gdiss_lim")

    @pytest.mark.usefixtures("tmodel")
    def test_phi_dict(self, tmodel):
        phi_dict = tmodel.phi_dict
        test_phi_dict = {'c': {'c': Q_(0, 'V'), 'm': Q_(-0.16, 'V'), 'e': Q_(0.06, 'V')},
                         'm': {'c': Q_(0.16, 'V'), 'm': Q_(0, 'V'), 'e': Q_(0, 'V')},
                         'e': {'c': Q_(-0.06, 'V'), 'm': Q_(0, 'V'), 'e': Q_(0, 'V')}}

        assert phi_dict == test_phi_dict

    @pytest.mark.usefixtures("tmodel")
    def test_rxn_thermo_params(self, tmodel):
        for rxn in tmodel.reactions:
            assert hasattr(rxn, "drG0")
            assert hasattr(rxn, "drG0prime")
            assert hasattr(rxn, "drGtransport")
            assert hasattr(rxn, "drGtransform")
            assert hasattr(rxn, "drG_h_transport")
            assert hasattr(rxn, "drG_c_transport")
            assert hasattr(rxn, "drG")
            assert hasattr(rxn, "drG_SE")

    @pytest.mark.usefixtures("tmodel")
    def test_met_thermo_params(self, tmodel):
        for met in tmodel.metabolites:
            assert hasattr(met, "upper_bound")
            assert hasattr(met, "lower_bound")
            assert hasattr(met, "concentration")
            assert hasattr(met, "accession")
            assert hasattr(met, "dfG0")
            assert hasattr(met, "dfG0prime")

    @pytest.mark.usefixtures("tmodel")
    def test_proton_dict(self, tmodel):
        proton_dict = tmodel.proton_dict
        assert len(proton_dict) == len(tmodel.compartments)

    @pytest.mark.usefixtures("tmodel")
    def test_charge_dict(self, tmodel):
        charge_dict = tmodel.charge_dict
        assert len(charge_dict) == len(tmodel.compartments)

    @pytest.mark.usefixtures("tmodel")
    def test_mg_dict(self, tmodel):
        mg_dict = tmodel.mg_dict
        assert len(mg_dict) == len(tmodel.compartments)

    #test compartment parents functionality
    @pytest.mark.usefixtures("tmodel")
    def test_compartment_parents(self, tmodel):
        #check default compartment parent keys match compartments
        for comp in tmodel.compartments:
            assert comp in tmodel.compartment_parents.keys()

        #test setting compartment parents
        tmodel.compartment_parents['m'] = 'c'
        assert tmodel.compartment_parents['m'] == 'c'
        #check that inner compartments are reset
        assert tmodel._inner_compartments is None
        #check that inner compartments are calculated correctly
        assert tmodel.inner_compartments[('c','m')] == 'm'


class TestExcelModel:

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
