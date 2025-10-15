#! 


# Packages needed:
import unittest
from thermo_flux.io import load_excel as ex
from thermo_flux.io import load_gams as gs
from thermo_flux.io import helper_load as hl
from thermo_flux.io import load_sbml as ip
from thermo_flux.core.model import ThermoModel, ThermoReaction, ThermoMetabolite
import numpy as np
import pandas as pd
import _pickle as pickle
import collections
import h5py
import cobra


class TestReadExcelYeast(unittest.TestCase):

    CASE = "yeast"
    PATH_FILE = "examples\\yeast\\yeast_v3.SV_ENS_03.xlsx"
    PATH_KEGG = "examples\\yeast\\yeast_kegg-inchi_id.csv"
    PATH_GAMS = "examples\\yeast\\model_yeast_from-GAMS.xlsx"
    EDIT_METS = {"ficytc": {"charge": 3},
                 "focytc": {"charge": 2},
                 "ac": {"charge": -1, "formula": "C2H3O2"},
                 "h2o": {"formula": "H2O"},
                 "succ": {"charge": -2}}

    @classmethod
    def setUpClass(cls):
        cls.model = ex.create_thermo_model(cls.CASE,
                                          model_excel=cls.PATH_FILE,
                                          keggids_csv=cls.PATH_KEGG,
                                          edit_mets=cls.EDIT_METS)

        cls.gams = hl.excel_to_df(cls.PATH_GAMS)

        cls.maxDiff = None

    def test_model_class(self):
        to_test = self.model
        self.assertTrue(isinstance(to_test, ThermoModel))
        return

    def test_rxn_class(self):
        to_test = self.model.reactions[0]
        self.assertTrue(isinstance(to_test, ThermoReaction))
        return

    def test_met_class(self):
        to_test = self.model.metabolites[0]
        self.assertTrue(isinstance(to_test, ThermoMetabolite))
        return

    def test_retrieve_metidx_single(self):
        idx_list = hl._retrieve_metabolite_idx_from_name(self.model, "ac")
        name_list = [self.model.metabolites[i].id.rsplit("_", 1)[0] for i in idx_list]
        no_comps = len(self.model.compartments)
        self.assertEqual(set(name_list), {"ac"})
        self.assertTrue(len(idx_list) <= no_comps * 1)
        return

    def test_retrieve_met_idx_list(self):
        idx_list = hl._retrieve_metabolite_idx_from_name(self.model, ["ac", "etoh"])
        name_list = [self.model.metabolites[i].id.rsplit("_", 1)[0] for i in idx_list]
        no_comps = len(self.model.compartments)
        self.assertEqual(set(name_list), {"ac", "etoh"})
        self.assertTrue(len(idx_list) <= no_comps * 2)
        return

    ############
    ### SETS ###
    ############
    def test_met_names(self):
        fromgams = set(self.gams["mets"].index)
        fromgams.add("biomass")
        fromgams.remove("electron")
        frommodel = set([m.id.rsplit("_", 1)[0] for m in self.model.metabolites])
        self.assertEqual(fromgams, frommodel)
        return

    def test_comps(self):
        fromgams = set(self.gams["comps"].index)
        frommodel = set(self.model.compartments)
        self.assertEqual(fromgams, frommodel)
        return 

    def test_rxns(self):
        fromgams = set(self.gams["rxns"].index.str.lower())
        fromgams.remove("biomass")
        frommodel = set([rxn.id.lower().replace("int", "") for rxn in self.model.reactions
                         if "biomass" not in rxn.id])
        self.assertEqual(fromgams, frommodel)
        return

    def test_extra_biomass_rxns(self):
        frommodel = set([rxn.id.lower() for rxn in self.model.reactions if "biomass" in rxn.id])
        biomass_rxns = {"biomass_c", "biomass_ce", "biomass_ex"}
        return self.assertEqual(biomass_rxns, frommodel)

    def test_extra_biomass_mets(self):
        frommodel = set([m.id for m in self.model.metabolites if "biomass" in m.id])
        biomass_mets = {"biomass_c", "biomass_e"}
        return self.assertEqual(biomass_mets, frommodel)

    def test_met_kegg(self):
        frommodel = set([m.annotation["kegg"] for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "ac"])
        return self.assertEqual({"C00033"}, frommodel)

    def test_met_kegg_biomass(self):
        frommodel = set([m.annotation["metanetx.chemical"] for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "biomass"])
        return self.assertEqual({"BIOMASS"}, frommodel)

    # def test_met_kegg_protein(self):
    #     frommodel = set([m.annotation["metanetx.cheumical"] for m in self.model.metabolites
    #                      if m.name.rsplit("_", 1)[0] == "protein"])
    #     return self.assertEqual({"UNKNOWN"}, frommodel)

    def test_met_kegg_charge(self):
        frommodel = set([m.annotation["metanetx.chemical"] for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "charge"])
        return self.assertEqual({"MNXM45842"}, frommodel)

    def test_met_inchi_ubiquinol(self):
        frommodel = set([m.annotation["InChI"] for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "q8h2"])
        return self.assertEqual({"InChI=1S/C49H78O4/c1-36(2)20-13-21-37(3)22-14-23-38(4)24-15-25-39(5)26-16-27-40(6)28-17-29-41(7)30-18-31-42(8)32-19-33-43(9)34-35-45-44(10)46(50)48(52-11)49(53-12)47(45)51/h20,22,24,26,28,30,32,34,46-47,50-51H,13-19,21,23,25,27,29,31,33,35H2,1-12H3"}, frommodel)

    def test_electron_removed(self):
        frommodel = set([m for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "electron"])
        return self.assertEqual(set(), frommodel)

    def test_int_transportrxns_removed(self):
        frommodel = set([rxn.id for rxn in self.model.reactions
                         if rxn.id.startswith("int")])
        return self.assertEqual(set(), frommodel)

    def test_proton_formula(self):
        frommodel = set([m.formula for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "h"])
        return self.assertEqual({"H1"}, frommodel)

    def test_charge_charge(self):
        frommodel = set([m.charge for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "charge"])
        return self.assertEqual({1}, frommodel)

    def test_special_mets_setcharge(self):
        frommodel = set([m.charge for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "succ"])
        return self.assertEqual({-2}, frommodel)

    def test_special_mets_setformula(self):
        frommodel = set([m.formula for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "h2o"])
        return self.assertEqual({"H2O"}, frommodel)

    def test_special_mets_settwoparams_formula(self):
        frommodel = set([m.formula for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "ac"])
        return self.assertEqual({"C2H3O2"}, frommodel)

    def test_special_mets_settwoparams_charge(self):
        frommodel = set([m.charge for m in self.model.metabolites
                         if m.name.rsplit("_", 1)[0] == "ac"])
        return self.assertEqual({-1}, frommodel)

    ##################
    ### Parameters ###
    ##################

    # # TO-DO: how to check this? not included in model_Ecoli.gdx
    # def test_pH(self):   # Parameters sheet
    #     return

    # def test_membrpot(self):
    #     return

    def test_temperature(self):   # Parameters sheet
        fromgams = self.gams["Scalar_param"].loc["Temp", "Value"]
        frommodel, = [val for k, val in self.model._parameters.items() if "temp" in k.lower()]
        self.assertEqual(fromgams, frommodel)
        return

    def test_gasconst(self):   # Parameters sheet
        fromgams = self.gams["Scalar_param"].loc["GasCons", "Value"]
        frommodel, = [val for k, val in self.model._parameters.items() if "gas" in k.lower()]
        self.assertEqual(fromgams, frommodel)
        return

    # def test_transmembrane(self):  # Transmembrane reactions sheet
    #     return

    def test_S_beforebalancing(self):
        fromgams = _build_Sdict_gams(self.gams)
        frommodel = _build_Sdict_cobra(self.model)

        bm_rxns_model = [rxn.id.lower() for rxn in self.model.reactions if "biomass" in rxn.id]
        bm_rxns_model += ["ex_c", "cex", "ionpump", "ionpumpm", "poscharge_ex", "posqt", "posqtm"]  
        frommodel = _pop_list_from_dict(frommodel, bm_rxns_model)

        bm_rxns_gams = ["biomass", "ex_c", "cex", "ionpump", "ionpumpm", "poscharge_ex", "posqt", "posqtm"]
        fromgams = _pop_list_from_dict(fromgams, bm_rxns_gams)

        return self.assertDictEqual(fromgams, frommodel)


class TestReadExcelEcoli(TestReadExcelYeast):
    CASE = "ecoli"
    PATH_FILE = "examples\\ecoli\\model.xlsx"
    PATH_KEGG = "examples\\ecoli\\ecoli_kegg-inchi_id.csv"
    PATH_GAMS = "examples\\ecoli\\model_Ecoli_from-gams.xlsx"
    EDIT_METS = {"ac": {"charge": -1, "formula": "C2H3O2"},
                 "h2o": {"formula": "H2O"},
                 "succ": {"charge": -2}}


class TestReadGamsYeast(unittest.TestCase):
    CASE = "yeast"
    base_path = base_path = "examples\\yeast\\"
    PATH_FILE = f"{base_path}model_yeast_from-GAMS.xlsx"
    PATH_KEGG = f"{base_path}yeast_kegg-inchi_id.csv"
    PATH_REGR = f"{base_path}yeast_regdata_from-GAMS.xlsx"
    PATH_VA = f"{base_path}yeast_vabounds_from-GAMS.xlsx"
    PATH_GAMS = f"{base_path}model_yeast_from-GAMS.xlsx"
    EDIT_METS = {"ac": {"charge": -1, "formula": "C2H3O2"},
                 "h2o": {"formula": "H2O"},
                 "succ": {"charge": -2}}
    RXN_EX_ACE = "ac_EX"
    RXN_EX_GLC = "glc-D_EX"

    @classmethod
    def setUpClass(cls):
        cls.model = gs.create_thermo_model_from_gams(cls.CASE,
                                                    model_excel=cls.PATH_FILE,
                                                    keggids_csv=cls.PATH_KEGG,
                                                    regr_excel=cls.PATH_REGR,
                                                    va_excel=cls.PATH_VA,
                                                    edit_mets=cls.EDIT_METS)
        cls.gams = hl.excel_to_df(cls.PATH_GAMS)
        cls.maxDiff = None

    def test_model_class(self):
        to_test = self.model
        self.assertTrue(isinstance(to_test, ThermoModel))
        return

    def test_rxn_class(self):
        to_test = self.model.reactions[0]
        self.assertTrue(isinstance(to_test, ThermoReaction))
        return

    def test_met_class(self):
        to_test = self.model.metabolites[0]
        self.assertTrue(isinstance(to_test, ThermoMetabolite))
        return

    def test_lnc_colnames(self):
        frommodel = set(self.model.lncLimits.columns)
        colnames = {"lnc_lo", "lnc_up", "va_lnc_lo", "va_lnc_up"}
        return self.assertEqual(colnames, frommodel)

    def test_flux_colnames(self):
        frommodel = set(self.model.fluxLimits.columns)
        colnames = {"v_lo", "v_up", "va_v_lo", "va_v_up"}
        return self.assertEqual(colnames, frommodel)

    def test_drg_colnames(self):
        frommodel = set(self.model._thermo.columns)
        colnames = {"drG0", "drGSE", "drGerror", "drG_lo", "drG_up",
                    "va_drG_lo", "va_drG_up"}
        return self.assertEqual(colnames, frommodel)

    def test_nullspace_dimensions(self):
        frommodel = self.model._nullspace
        rows, cols = frommodel.shape

        fromgams = self.gams["K"]
        rows_gams = len(fromgams.index.unique())
        cols_gams = len(fromgams.dim2.unique())

        self.assertEqual(rows, rows_gams)
        self.assertEqual(cols, cols_gams)
        return

    # TO-DO: in the case of E.coli, the coefficients are not all integers
    def test_nullspace_integer(self):
        frommodel = self.model._nullspace
        k_as_column = frommodel.melt().dropna()["value"]
        comp = k_as_column.apply(lambda x: float(x).is_integer())
        return self.assertTrue(all(comp))

    def test_find_met_exchanged(self):
        m = self.model
        self.assertEqual("biomass_e", gs._find_met_exchanged(m, "biomass_EX"))
        self.assertEqual("ac_e", gs._find_met_exchanged(m, self.RXN_EX_ACE))
        self.assertEqual("glc-D_e", gs._find_met_exchanged(m, self.RXN_EX_GLC))
        self.assertRaises(ValueError, gs._find_met_exchanged, m, "PYK")
        return

    # TO-DO: this doesn't work because some of the bounds in the input files are actually not narrower 
    # than the defaults (e.g. pyr_c, pi_c, ...)
    # def test_VAdrgbounds_narrower(self):
    #     frommodel = self.model._thermo[["v_lo", "va_v_lo", "v_up", "va_v_up"]].dropna()
    #     range_default = frommodel["v_up"] - frommodel["v_lo"]
    #     range_va = frommodel["va_v_up"] - frommodel["va_v_lo"]
    #     return self.assertTrue(all(range_va <= range_default))

    # def test_VAlncbounds_narrower(self):
    #     frommodel = self.model.lncLimits[["v_lo", "va_v_lo", "v_up", "va_v_up"]].dropna()
    #     range_default = frommodel["v_up"] - frommodel["v_lo"]
    #     range_va = frommodel["va_v_up"] - frommodel["va_v_lo"]
    #     return self.assertTrue(all(range_va <= range_default))

    # def test_VAdrgbounds_narrower_up(self):
    #     frommodel = self.model._thermo[["v_up", "va_v_up"]].dropna()
    #     comp = frommodel["va_v_up"] <= frommodel["v_up"]
    #     return self.assertTrue(all(comp))

    # def test_VAdrgbounds_narrower_lo(self):
    #     frommodel = self.model._thermo[["v_lo", "va_v_lo"]].dropna()
    #     comp = frommodel["va_v_lo"] >= frommodel["v_lo"]
    #     return self.assertTrue(all(comp))

    # def test_VAlncbounds_narrower_up(self):
    #     frommodel = self.model.lncLimits[["v_up", "va_v_up"]].dropna()
    #     comp = frommodel["va_v_up"] <= frommodel["v_up"]
    #     return self.assertTrue(all(comp))

    # def test_VAlncbounds_narrower_lo(self):
    #     frommodel = self.model.lncLimits[["v_lo", "va_v_lo"]].dropna()
    #     comp = frommodel["va_v_lo"] >= frommodel["v_lo"]
    #     return self.assertTrue(all(comp))


class TestReadGamsEcoli(TestReadGamsYeast):
    CASE = "ecoli"
    base_path = "examples\\ecoli\\"
    PATH_FILE = f"{base_path}model_Ecoli_from-gams.xlsx"
    PATH_KEGG = f"{base_path}ecoli_kegg-inchi_id.csv"
    PATH_REGR = f"{base_path}ecoli_regdata_from-gams.xlsx"
    PATH_VA = f"{base_path}ecoli_vabounds_from-gams.xlsx"
    PATH_GAMS = f"{base_path}model_Ecoli_from-gams.xlsx"
    EDIT_METS = {"ficytc": {"charge": 3},
                 "focytc": {"charge": 2},
                 "ac": {"charge": -1, "formula": "C2H3O2"},
                 "h2o": {"formula": "H2O"},
                 "succ": {"charge": -2}}
    RXN_EX_ACE = "EX_ac"
    RXN_EX_GLC = "EX_glc"


class TestReadSBML(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open("thermo_flux\\io\\testpickled_yeast-from-GAMS.pkl", "rb") as input:
            model_original = pickle.load(input)
        cls.model = model_original

        cls.imported = ip.read_from_sbml("thermo_flux\\io\\testsbml_yeast-from-GAMS.sbml")
        cls.maxDiff = None

    def test_model_class(self):
        to_test = self.imported
        self.assertTrue(isinstance(to_test, ThermoModel))
        return

    def test_rxn_class(self):
        to_test = self.imported.reactions[0]
        self.assertTrue(isinstance(to_test, ThermoReaction))
        return

    def test_met_class(self):
        to_test = self.imported.metabolites[0]
        self.assertTrue(isinstance(to_test, ThermoMetabolite))
        return

    #############
    ### Utils ###
    #############
    def test_dictconverter_allfloats(self):
        dict_in = {"a": "1", "b": "2.5"}
        dict_out = {"a": 1.0, "b": 2.5}
        dict_in_conv = ip._convert_entries_to_float(dict_in)
        self.assertDictEqual(dict_in_conv, dict_out)
        return

    def test_dictconverter_stringsfloats(self):
        dict_in = {"a": "1", "b": "B"}
        dict_out = {"a": 1.0, "b": "B"}
        dict_in_conv = ip._convert_entries_to_float(dict_in)
        self.assertDictEqual(dict_in_conv, dict_out)
        return

    def test_checkisiterable_true(self):
        obj_in = [0, 1, 2]
        self.assertTrue(ip._check_is_iterable(obj_in))
        return

    def test_checkisiterable_false(self):
        obj_in = None
        self.assertFalse(ip._check_is_iterable(obj_in))
        return

    def test_noteswrapper_model(self):
        notes_imp = self.imported.notes
        notes_orig = self.model._parameters
        self.assertDictEqual(notes_imp, notes_orig)
        return

    def test_noteswrapper_reactions(self):
        idx = 0
        rxn_id = self.imported.reactions[idx].id
        notes_imp = self.imported.reactions[idx].notes
        notes_orig = self.model.fluxLimits.loc[rxn_id].to_dict()
        notes_orig_2 = self.model._thermo.loc[rxn_id].to_dict()
        notes_orig.update(notes_orig_2)
        # self.assertDictEqual(notes_imp, notes_orig)
        np.testing.assert_equal(notes_imp, notes_orig)
        return

    def test_noteswrapper_metabolites(self):
        idx = 0
        met_id = self.imported.metabolites[idx].id
        notes_imp = self.imported.metabolites[idx].notes
        notes_orig = self.model.lncLimits.loc[met_id].to_dict()
        # self.assertDictEqual(notes_imp, notes_orig)
        np.testing.assert_equal(notes_imp, notes_orig)
        return

    def test_notes_reactions_drg0(self):
        fromsbml = self.imported
        value_sbml = _build_dict_notes(fromsbml.reactions, "drG0")
        original = self.model._thermo["drG0"].to_dict()
        original = {k: val for k, val in original.items() if "biomass" not in k and "poscharge" not in k}
        np.testing.assert_equal(value_sbml, original)
        return

    def test_notes_reactions_drgSE(self):
        fromsbml = self.imported
        value_sbml = _build_dict_notes(fromsbml.reactions, "drGSE")
        original = self.model._thermo["drGSE"].to_dict()
        original = {k: val for k, val in original.items() if "biomass" not in k and "poscharge" not in k}
        np.testing.assert_equal(value_sbml, original)
        return

    def test_notes_reactions_drgerror(self):
        fromsbml = self.imported
        value_sbml = _build_dict_notes(fromsbml.reactions, "drGerror")
        original = self.model._thermo["drGerror"].to_dict()
        original = {k: val for k, val in original.items() if "biomass" not in k and "poscharge" not in k}
        np.testing.assert_equal(value_sbml, original)
        return

    def test_notes_metabolite_lncup(self):
        fromsbml = self.imported
        value_sbml = _build_dict_notes(fromsbml.metabolites, "lnc_up")
        original = self.model.lncLimits["lnc_up"].to_dict()
        np.testing.assert_equal(value_sbml, original)
        return

    def test_notes_metabolite_lncvaup(self):
        fromsbml = self.imported
        value_sbml = _build_dict_notes(fromsbml.metabolites, "va_lnc_up")
        original = self.model.lncLimits["va_lnc_up"].to_dict()
        np.testing.assert_equal(value_sbml, original)
        return


# # class TestReadHDF5(unittest.TestCase):
# #     @classmethod
# #     def setUpClass(cls):

# #         with open("thermo_flux\\io\\test_model_pickled.pkl", "rb") as input:
# #             model_original = pickle.load(input)
# #         cls.model = model_original

# #         test_model = cobra.Model("ecoli")
# #         test_model._thermo = cls.model._thermo.copy()
# #         test_model._thermo = 0   # mock case for conflicting values
# #         cls.test_model = test_model.copy()
# #         cls.filename = "thermo_flux\\io\\test_model.hdf5"
# #         cls.imported = ip.read_from_hdf5(test_model, cls.filename, priority="hdf5")

# #         cls.maxDiff = None

# #     # HIGH-LEVEL FUNCTIONS
# #     def test_fetch_df_type(self):
# #         from_import = self.imported._stoichmatrix
# #         original = self.model._stoichmatrix
# #         self.assertEqual(type(from_import), type(original))
# #         return

# #     def test_fetch_df_shape(self):
# #         from_import = self.imported._stoichmatrix
# #         original = self.model._stoichmatrix
# #         self.assertEqual(from_import.shape, original.shape)
# #         return

# #     def test_fetch_df_values(self):
# #         from_import = self.imported._stoichmatrix
# #         original = self.model._stoichmatrix
# #         self.assertTrue(np.allclose(from_import.values, original.values, equal_nan=True))
# #         return

# #     def test_fetch_params(self):
# #         from_import = self.imported._parameters
# #         original = self.model._parameters
# #         self.assertDictEqual(from_import, original)
# #         return

# #     # LOWER-LEVEL FUNCTIONS
# #     # def test_solveconflicts_missing_attribute(self):
# #     #     with h5py.File(self.filename, "r") as f:
# #     #         f_ = f[self.test_model.id]
# #     #         data_grp = f_["data"]
# #     #         args = (self.test_model, data_grp, "hdf5")
# #     #         self.assertRaises(AttributeError, ip._solve_conflicts, *args)
# #     #     return

# #     def test_solveconflicts_priority_hdf5(self):
# #         from_import = self.imported._thermo
# #         original = self.model._thermo
# #         self.assertTrue(np.allclose(from_import.values, original.values, equal_nan=True))
#         return


#############
### UTILS ###
#############

def _build_Sdict_cobra(cobramodel):
    sdict = {}
    for r in cobramodel.reactions:
        rname = r.id.lower()  # .replace("int", "")

        try:
            sdict[rname]
        except KeyError:
            sdict[rname] = {}

        stoich_rxn = {met.id: float(s) for met, s in r.metabolites.items() if met.name.rsplit("_", 1)[0] not in ["h", "charge"]}

        counter = collections.Counter()   # https://www.geeksforgeeks.org/python-sum-list-of-dictionaries-with-same-key/
        for d in [sdict[rname], stoich_rxn]: 
            counter.update(d)

        sdict[rname] = dict(counter)

        try:
            mets_to_drop = [met for met, s in sdict[rname].items() if s==0]
            sdict[rname] = _pop_list_from_dict(sdict[rname], mets_to_drop)
        except TypeError:
            pass

    return sdict

# TO-DO: better reconstruction of the S-matrix, including charge-transport reactions,
# while still ignoring charge balance in the remaining
def _build_Sdict_gams(df_gams):
    df = df_gams["S"].reset_index()
    df["dim3"] = df["dim3"].str.lower()
    df["met"] = df["dim1"] + "_" + df["dim2"]    # metaboliteName_compartment
    sdict = df.groupby("dim3").apply(
                                    lambda x: {row["met"]: float(row["Value"]) for _, row in x.iterrows() 
                                    if row["met"].rsplit("_", 1)[0] not in ["h", "charge"]}  # this reconstructs the S-matric but ignores all charges (including in Cex and EX_C!)
                                    ).to_dict()
    return sdict


def _build_dict_notes(iterator, key):
    dict_out = {}
    for x in iterator:
        if key in x.notes.keys():
            dict_out[x.id] = x.notes[key]
        else:
            pass
    return dict_out


def _pop_list_from_dict(inputdict, poplist):
    for p in poplist:
        try:
            inputdict.pop(p)
        except KeyError:
            pass
    return inputdict
