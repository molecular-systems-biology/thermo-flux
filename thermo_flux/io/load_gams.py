#! 


# Packages needed:
from cobra import Model, Reaction, Metabolite
# from . import helper_load as hl
from thermo_flux.io import helper_load as hl
from thermo_flux.core.model import ThermoModel
import numpy as np

def create_S_matrix_from_gams(modelname, model_excel, keggids_csv,
                                  regr_excel=None, va_excel=None,
                                  edit_mets=None, kwargs_rxns={}, **kwargs):
    """
    Creates and returns a ThermoModel, populated with the information
    contained in an excel file, generated from a GAMS model, .gdx.
    """
    model = Model(modelname)

    df_map = read_inputs(model_excel)
    df_map_regr = read_inputs(regr_excel)
    df_map_va = read_inputs(va_excel)

    model = add_reactions_from_Smatrix(model, df_map, **kwargs_rxns)

    thermo_model = ThermoModel(model, **kwargs)

    return thermo_model

def create_thermo_model_from_gams(modelname, model_excel, keggids_csv,
                                  regr_excel=None, va_excel=None,
                                  edit_mets=None, kwargs_rxns={}, **kwargs):
    """
    Creates and returns a ThermoModel, populated with the information
    contained in an excel file, generated from a GAMS model, .gdx.
    """
    model = Model(modelname)

    df_map = read_inputs(model_excel)
    df_map_regr = read_inputs(regr_excel)
    df_map_va = read_inputs(va_excel)

    model = add_reactions_from_Smatrix(model, df_map, **kwargs_rxns)

   # model = hl.add_biomass_reactions(model)

    model = update_metabolite_info(model, keggids_csv)

    model = hl.edit_special_metabolites(model, edit_mets)

    model = add_thermodynamic_data(model, df_map, df_map_regr, df_map_va)

    model = add_concentration_data(model, df_map, df_map_va)

    model = add_flux_data(model, df_map, df_map_va)

    model = add_nullspace(model, df_map)

    model = add_scalar_params(model, df_map)

    thermo_model = ThermoModel(model, **kwargs)

    return thermo_model


def read_inputs(inputpath):
    """
    Returns a dictionary mapping the name of each sheet on the input
    excel file to the dataframe containing the data, or None, if the
    path to the file is not specified.
    """
    if inputpath is not None:
        df_map = hl.excel_to_df(inputpath)
    else:
        print("input is none!")
        df_map = None
    return df_map


def add_reactions_from_Smatrix(model, df, **kwargs):
    """
    Takes the input dataframe containing the data from GAMS and populates
    the COBRA model with the reactions and metabolites specified in the
    stoichiometry matrix.
    """
    print("*** Building model from stoichiometrix matrix ***")

    df_S = _build_stoich_matrix(df)

    dict_mets = _get_met_dictionary(df_S)

    model = _add_reactions_to_model(model, df_S, dict_mets, **kwargs)

    model._stoichmatrix = df_S

    return model


def _build_stoich_matrix(df):
    """
    Returns a dataframe where the row are indexed by the reaction names,
    and the columns are the compartment-specific metabolites.
    The values of the matrix entries are the stoichiometric coefficients
    of each metabolite in the given reaction.

    S = [rxns x mets]
    """
    df_S = df["S"].reset_index()
    df_S.rename(columns={"dim1": "mets",
                         "dim2": "comps",
                         "dim3": "rxns",
                         "Value": "coeff"}, inplace=True)
    df_S["mets"] = df_S["mets"] + "_" + df_S["comps"]
    df_S = df_S.pivot(columns="mets", values="coeff", index="rxns")
    return df_S


def _get_met_dictionary(df):
    """
    Returns a dictionary of metabolite_id: cobra.Metabolite()
    """
    dict_mets = {}
    for met_id in df.columns:
        met_name, met_comp = met_id.rsplit("_", 1)
        dict_mets[met_id] = Metabolite(met_id,
                                       name=met_name,
                                       compartment=met_comp)
    return dict_mets


def _add_reactions_to_model(model, df_S, dict_mets, lo=-500, up=500):
    """
    Adds the reactions to the model according to the stoichiometry defined
    in df_S.
    Sets the upper and lower bounds for the fluxes, by default +/- 500.
    """
    for rxn_id, row in df_S.iterrows():
        reaction = Reaction(rxn_id,
                            name=rxn_id,
                            lower_bound=lo,
                            upper_bound=up)

        rxn_stoich = {}
        row_dict = row.dropna().to_dict()
        for met_id, coeff in row_dict.items():
            met = dict_mets[met_id]
            rxn_stoich[met] = coeff

        reaction.add_metabolites(rxn_stoich)

        model.add_reactions([reaction])

    return model


def update_metabolite_info(model, keggids_csv):
    """
    Updates metabolite information: name, abbreviation, compartment and KEGG ID.
    """
    print("*** Updating metabolite information ***")
    kegg_ids = hl.read_kegg_csv(keggids_csv)
    for met in model.metabolites:
        met.charge = 0
        met = hl.add_kegg(met, met.name, kegg_ids)
    return model


def add_thermodynamic_data(model, df, df_regr, df_va):

    model = _add_drg0(model, df)
    model = _update_drg0_exchanges(model, df)
    model = _update_drg0_biosyntesis(model, df, df_regr)
    model = _add_drgerror(model, df_regr)
    model = _add_limits_wrapper(_add_drglimits, model, df, df_va, df_regr=None)
    return model


# TO-DO: use df["rxns"] or model.reactions to set the index? (check new biomass reactions)
def _add_drg0(model, df):
    """
    Adds an attribute <thermo> to the model, containing the value of drG0
    and its uncertainty for each reaction.
    """
    df_drg0 = df["drGt0tr"].reindex(df["rxns"].index, fill_value=0)
    df_drg0.rename(columns={"Value": "drG0"}, inplace=True)
    df_drg0["drGSE"] = df["drGSE"]["Value"]
    df_drg0.drop([x for x in df_drg0.columns if "Unnamed" in x],
                 axis=1, inplace=True)

    model._thermo = df_drg0

    return model


def _update_drg0_exchanges(model, df):
    """
    Update the values of drG0 for the exchange reactions with the dfG0
    of the metabolite being exchanged.
    """
    ex_rxns = [r.id for r in model.reactions if "_ex" in r.id.lower() or "ex_" in r.id.lower()]

    for rxn in ex_rxns:
        met_id = _find_met_exchanged(model, rxn)
        m, c = met_id.rsplit("_", 1)

        df_dfg0 = df["dfG0"].reset_index()
        try:
            dfg, = df_dfg0.loc[(df_dfg0["dim1"] == m) & (df_dfg0["dim2"] == c), "Value"]
        except ValueError:
            dfg = 0

        model._thermo.loc[rxn, "drG0"] = dfg

    return model


# TO-DO: worth doing this when the model does not 
def _update_drg0_biosyntesis(model, df, df_regr):
    """
    For models in which the biomass production is subdivided in the
    production of its constituent macromolecules, drG0 must be calculated
    for each of this biosynthetic reactions.
    """
    try:
        df_bio = df_regr["dfG_bmm"]
    except KeyError as e:
        print(f"Provided regression data does not have a sheet called {e}")
        df_bio = None
    except TypeError:
        print("Did not provide regression data")
        df_bio = None

    if df_bio is not None:
        biosynth_rxns = df_bio.loc[df_bio["Value"] != 0]  # .index.to_list()

        for rxn, row_bmm in biosynth_rxns.iterrows():
            stoich = model._stoichmatrix.loc[rxn].dropna()
            met_ids = stoich.reset_index()["mets"]

            dfg0_all = df["dfG0"].reset_index()
            dfg0_all["mets"] = dfg0_all["dim1"] + "_" + dfg0_all["dim2"]

            dfg0_mets = dfg0_all.loc[dfg0_all.mets.isin(met_ids)].set_index("mets")
            dfg0_macromol = row_bmm["Value"]

            drg0_bio = (dfg0_mets["Value"] * stoich).sum() + dfg0_macromol

            model._thermo.loc[rxn, "drG0"] = drg0_bio

    return model


def _add_drgerror(model, df_regr):
    """
    Reads the values of drGerror from the input excel file, if provided.
    Values stored in the <thermo> attribute of cobra.Model.
    """
    try:
        model._thermo["drGerror"] = df_regr["drGerror"]["Value"]
    except KeyError as e:
        print(f"Provided regression data does not contain {e}...")
    except TypeError:
        print("Did not provide regression data as input...")
    finally:
        model._thermo["drGerror"] = np.nan
        print("... no drGerror values added.")

    return model


def _add_limits_wrapper(func, model, df, df_va, df_regr):
    """
    Wrapper for adding drGLimits to the model. Handles the different
    kinds of inputs.
    """
    model = func(model, df)

    if df_va is not None:
        model = func(model, df_va, prefix="va_")

    if df_regr is not None:
        model = func(model, df_regr, prefix="regr_")

    return model


def _add_drglimits(model, df, prefix=""):
    """
    Read and add the values of drGLimits (or dfGLimits, in the case of
    exchange reactions) to the dataframe storing the thermodynamic information.
    """
    try:
        drgLim = df["drGLimits"].pivot(columns="dim2", values="Value").fillna(0)
        drgLim.rename(columns={"lo": "drG_lo", "up": "drG_up"}, inplace=True)
        model._thermo[prefix+"drG_lo"] = drgLim["drG_lo"]
        model._thermo[prefix+"drG_up"] = drgLim["drG_up"]

        drgLim_ex = df["dfGLimits"].pivot(columns="dim2", values="Value").fillna(0)
        drgLim_ex.rename(columns={"lo": "drG_lo", "up": "drG_up"}, inplace=True)
        model._thermo[prefix+"drG_lo"].update(drgLim_ex["drG_lo"])
        model._thermo[prefix+"drG_up"].update(drgLim_ex["drG_up"])

    except KeyError as e:
        print(f"Provided V.A. data does not contain {e}...")
    except TypeError:
        print("Did not provide V.A. data as input...")

    return model


def _find_met_exchanged(model, rxn):
    """
    Returns the ID of a metabolite that is part of an exchange reaction
    """
    stoich = model.reactions.get_by_id(rxn).metabolites
    met_id, = [met.id for met, _ in stoich.items()]
    return met_id


def add_concentration_data(model, df, df_va):
    """
    Adds new attribute to the cobra model, containing a dataframe that
    describes the upper and lower bounds of metabolite concentration (lnc).
    Metabolites are identified by name_compartment.
    """
    model = _add_limits_wrapper(_add_lnclimits, model, df, df_va, df_regr=None)

    return model


def _add_lnclimits(model, df, prefix=""):
    """
    Read and add the lncLimits to an attribute of cobra.Model named <lncLimits>.
    """
    df_lnc = df["lncLimits"].reset_index()
    df_lnc.rename(columns={"dim1": "mets",
                           "dim2": "comps",
                           "dim3": "bndID",
                           }, inplace=True)
    df_lnc["mets"] = df_lnc["mets"] + "_" + df_lnc["comps"]
    df_lnc = df_lnc.pivot(columns="bndID", values="Value", index="mets")
    df_lnc.rename(columns={"lo": "lnc_lo", "up": "lnc_up"}, inplace=True)

    if hasattr(model, "lncLimits"):
        model.lncLimits[prefix+"lnc_lo"] = df_lnc["lnc_lo"]
        model.lncLimits[prefix+"lnc_up"] = df_lnc["lnc_up"]
    else:
        model.lncLimits = df_lnc

    return model


def add_flux_data(model, df, df_va):
    """
    Adds new attribute to the cobra model, containing a dataframe that
    describes the upper and lower bounds of fluxes of each reaction.
    """
    model = _add_limits_wrapper(_add_fluxlimits, model, df, df_va, df_regr=None)
    return model


def _add_fluxlimits(model, df, prefix=""):
    """
    Read and add the FluxLimits to an attribute of cobra.Model named <fluxLimits>.
    """
    try:
        vlim = df["FluxLimits"].pivot(columns="dim2", values="Value").fillna(0)
        vlim.rename(columns={"lo": "v_lo", "up": "v_up"}, inplace=True)

        if hasattr(model, "fluxLimits"):
            model.fluxLimits[prefix+"v_lo"] = vlim["v_lo"]
            model.fluxLimits[prefix+"v_up"] = vlim["v_up"]
        else:
            model.fluxLimits = vlim

    except KeyError as e:
        print(f"Provided V.A. data does not contain {e}...")
    except TypeError:
        print("Did not provide V.A. data as input...")

    return model


def add_nullspace(model, df):
    """
    Adds new attribute to the cobra model, containing a dataframe with
    the nullspace vectors.
    """
    df_K = df["K"].reset_index()
    df_K.rename(columns={"dim1": "rxns",
                         "dim2": "nsID",
                         "Value": "coeff"}, inplace=True)
    df_K = df_K.pivot(columns="nsID", values="coeff", index="rxns")
    model._nullspace = df_K

    return model


def add_scalar_params(model, df):
    """
    Stores the scalar parameters in a dictionary, as an additional attribute
    to the cobra.Model.
    """
    df_params = df["Scalar_param"]
    params_dict = {g: row["Value"] for g, row in df_params.iterrows()}

    model._parameters = params_dict

    return model
