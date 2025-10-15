#! 


# Packages needed:
import re
import pandas as pd
import numpy as np
from cobra import Model, Reaction
# from . import helper_load as hl
from thermo_flux.io import helper_load as hl
from thermo_flux.core.model import ThermoModel


def create_thermo_model(modelname, model_excel, keggids_csv, edit_mets=None, **kwargs):
    """
    Creates and returns a ThermoModel, populated with the information
    contained in an excel sheet.
    The excel sheet should follow the template used by the MSB group.
    """
    model = Model(modelname)

    sheet_to_df_map = read_excelfile(model_excel)

    sheets_to_parse = ["Reactions", "Biomass Composition", "Exchange reactions",
                       "Transmembrane reactions", "Parameters", "Metabolites", ]
    assert set(sheets_to_parse).issubset(set(sheet_to_df_map.keys()))

    model = parse_sheet_data(model, sheet_to_df_map, sheets_to_parse)

    model = merge_transport_reactions(model)

    #model = hl.add_biomass_reactions(model)

    model = update_metabolite_info(model, keggids_csv)

    model = hl.edit_special_metabolites(model, edit_mets)

    model = delete_electrons(model)

    thermo_model = ThermoModel(model, **kwargs)

    return thermo_model


def read_excelfile(model_excel):
    """
    Returns a dictionary mapping the information contained in each of
    the sheets (of the excel file) to the sheet's name.
    """
    xls = pd.ExcelFile(model_excel)

    # List all sheets in the excel file
    print(xls.sheet_names)

    # Read all sheets to a map
    sheet_to_df_map = {}
    for sheet_name in xls.sheet_names:
        sheet_to_df_map[sheet_name] = xls.parse(sheet_name)

    return sheet_to_df_map


def parse_sheet_data(model, sheet_map, sheets_to_parse):
    """
    Retrieves the information from the input excel file, adding the
    information to a COBRA model.
    """
    for sheet_name in sheets_to_parse:
        print(f"*** Reading data from {sheet_name} ***")
        df = sheet_map[sheet_name]

        if sheet_name == "Parameters":
            model = _parse_parameters(model, df)

        elif sheet_name == "Transmembrane reactions":
            for _, row in df.iterrows():
                transmembrane, intracellular = _parse_transport_reactions(row)

                model = _add_reactions_to_model(model, row,
                                                rxn_str=transmembrane)  # add transmembrane part, default abbreviation

                if type(intracellular) is str:
                    int_rid = "int" + row["Abbrevation"]
                    model = _add_reactions_to_model(model, row,
                                                    rxn_str=intracellular,
                                                    rxn_id=int_rid)   # add intracellular part, modified rxn abbreviation

        elif sheet_name == "Metabolites":   # parsing the Metabolites sheet should come after defining all the model reactions (metabolites are defined there)
            for _, row in df.iterrows():
                model = _add_met_charge_and_formula(model, row)

        else:
            for _, row in df.iterrows():
                model = _add_reactions_to_model(model, row)

    return model


def _add_reactions_to_model(model, row, rxn_str=None, rxn_id=None):
    """
    Adds a reaction to the COBRA model, starting from a string which
    specifies the stoichiometry of the reaction.
    Mandatory argument: row of the excel sheet containing the information
    for the given reaction.
    Optionally takes as arguments strings explicitly specifying the
    reaction stoichiometry and ID.
    """
    if rxn_str is None:
        rxn_str = row["Reaction"]
    if rxn_id is None:
        rxn_id = row["Abbrevation"]   # note misspelling

    rxn = Reaction(rxn_id)
    model.add_reactions([rxn])
    rxn.build_reaction_from_string(rxn_str)

    return model


# TO-DO: remove hard-coded column names storing abbreviations, charge and formula
def _add_met_charge_and_formula(model, row):
    """ 
    Reads the charge and chemical formula from the input .xlsx
    file into the appropriate attributes of each Metabolite in the model.
    """

    met_abbreviation = row[0]
    
    try:
        met_charge = row["Unnamed: 8"]
    except KeyError:
        met_charge = None

    try:
        met_formula = row["Unnamed: 12"]
    except KeyError:
        met_formula = row["Formula (neutral)"]

    # if met_charge == np.nan:
    #     met_charge = None
    # if met_formula == np.nan:
    #     met_formula = None

    all_mets_to_edit = [met for met in model.metabolites if met.id[:-3] == str(met_abbreviation)]

    for met_to_edit in all_mets_to_edit:
        met_to_edit.charge = met_charge
        met_to_edit.formula = met_formula

    return model

# TO-DO: make more robust, handle different Exceptions
def _parse_transport_reactions(row):
    """
    Return strings with the stoichiometries of the transport reations,
    split into (1) transmembrane and (2) intracellular parts.
    These components are already discriminated in the original excel
    sheet, and are read from separate columns.
    """
    try:
        transmembrane = re.sub("{[^}]+}", "", row["Transmembrane transport"])
        intracellular = row["Intracellular part"]

    except KeyError:
        transmembrane = None
        intracellular = None

    return transmembrane, intracellular


def _parse_parameters(model, df):
    """
    Parses the data from the "Parameters" sheet in the excel file.
    Stores on the COBRA model a dictionary where the primary keys are
    the parameter names; the secondary keys, if they exist, specify the
    domains of validity for the values (e.g. the cell compartment).
    """
    params_dict = {}
    cols = None
    par = None

    for _, row in df.iterrows():

        rr = row.replace(np.nan, False).values
        rr_filt = rr[rr != False]
        indices = [i for i, val in enumerate(rr) if val in rr_filt]

        if len(indices) > 0:
            if 0 not in indices:
                cols = {k: val for k, val in zip(rr_filt, indices)}
            else:
                par = rr_filt[0]
        else:
            cols = {k: val for k, val in zip(rr_filt, indices)}
            par = None

        if par is not None:
            if len(cols) > 0:
                params_dict[par] = {key: rr[i] for key, i in cols.items()}
            else:
                params_dict[par] = rr[1]

    model._parameters = params_dict

    return model


def merge_transport_reactions(model):
    """
    Merges the intracellular and transmembrane parts of transport reactions.
    """
    to_merge = [rxn.id.strip("int") for rxn in model.reactions
                if rxn.id.startswith("int")]

    for rr in to_merge:
        id_int = "int" + rr
        id_trans = rr

        # Stoichiometry of the internal reaction
        stoich_int = model.reactions.get_by_id(id_int).metabolites

        # Update overall reaction stoichiometry
        model.reactions.get_by_id(id_trans).add_metabolites(stoich_int)

        # Pop internal reaction
        model.reactions.get_by_id(id_int).delete()

    return model


def update_metabolite_info(model, keggids_csv):
    """
    Updates metabolite information: name, abbreviation, compartment
    and KEGG ID.
    """
    print("*** Updating metabolite information ***")
    kegg_ids = hl.read_kegg_csv(keggids_csv)
    for met in model.metabolites:
        try:
            comp = re.search("(?<=\[)(.*)(?=\])", met.id).group()
            fullname = met.id
            basename = fullname[0:fullname.index("[")]
        except AttributeError:
            basename, comp = met.id.rsplit("_", 1)   # to handle biomass_c and biomass_e

        met.compartment = comp
        ID = basename + "_" + comp
        met.id = ID
        met.name = ID
        # met.charge = 0
        met = hl.add_kegg(met, basename, kegg_ids)

    return model


def delete_electrons(model):
    """
    Remove the electrons from the model.
    """
    print("*** Removing electrons from the model ***")
    e_idx = hl._retrieve_metabolite_idx_from_name(model, "electron")
    electrons = [model.metabolites[i] for i in e_idx]
    model.remove_metabolites(electrons)
    return model
