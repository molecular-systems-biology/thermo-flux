#! 


# Packages needed:
import cobra
import h5py
import itertools
import numpy as np
import pandas as pd
from thermo_flux.core.model import ThermoModel
from equilibrator_api import Q_


def load_model(sbml_filename, hdf5_filename=None, priority="sbml", **kwargs):
    """"

    """

    model = read_from_sbml(sbml_filename)

    if hdf5_filename is not None:
        model = read_from_hdf5(model, hdf5_filename, priority)

    thermo_model = ThermoModel(model, **kwargs)
    
    thermo_model = format_notes(thermo_model)

    return thermo_model


def read_from_hdf5(model, filename, priority):
    """
    Populates a COBRA model with data stored in a .hdf5 file, throwing
    warnings if there are conflicts with data already in the model.
    """

    with h5py.File(filename, "r") as f:

        f_ = f[model.id]

        # # Populate model parameters
        # params_grp = {"_parameters": f_["_parameters"]}
        # model = _solve_conflicts(model, params_grp, priority)

        # Populate data:
        data_grp = f_["data"]
        model = _solve_conflicts(model, data_grp, priority)

    return model


# TO-DO: leave last Exception?
def _solve_conflicts(model, from_hdf5, priority):
    """
    Compares the values stored in the instance of COBRA.Model with the
    values stored in the .hdf5 file.
    In the case of conflicting values for the same attribute, sets
    according to the user's preference, specified in the 'priority' argument.
    """

    for k in from_hdf5.keys():
        try:
            data_grp = from_hdf5[k]
            values = _fetch_input_data(data_grp, k)
            from_model = getattr(model, k)

            np.testing.assert_equal(from_model, values)

        except AttributeError:
            print(f"Attribute '{k}' not in the model. Loading values from .hdf5 file.")
            setattr(model, "_"+k, values)

        except AssertionError:
            print(f"Conflicting values for attribute '{k}'")
            print(f"Will give priority to the values from the .{priority} file.")

            if priority == "sbml":
                pass
            elif priority == "hdf5":
                setattr(model, "_"+k, values)

        except Exception:
            print(f"UNFORESEEN ERROR when reading '{k}' from the .hdf5 file")

    return model


def _fetch_input_data(group, key):
    """
    Retrieves data stored in a given group of the HDF5 file (name given
    by 'key'). When the attributes so indicate, the data is converted to
    a pandas.DataFrame.
    """
    attributes = dict(group.attrs.items())
    values = group.values()

    if len(values) == 0:
        values_out = attributes
    else:
        # if {"columns", "index"}.issubset(set(attributes.keys())):
        #     values_out = pd.DataFrame(group[key][:],
        #                               index=attributes["index"],
        #                               columns=attributes["columns"])

        if "units" in attributes.keys():
            values_out = Q_(group[key][:], units=group.attrs["units"])

        else:
            values_out = group[key][:]

    return values_out


def read_from_sbml(filename):
    """
    Reads an .sbml file and builds a COBRA model from it, ensuring that
    the 'notes' attributes of all relevant objects (Model, Reaction and
    Metabolite) have been converted to floating numbers, wherever possible.
    """

    model = cobra.io.read_sbml_model(filename)

    # model = format_notes(model)

    return model


def format_notes(model):
    """
    Formats the dictionaries stored in 'notes' of different Cobra objects
    (Model, Reaction and Metabolite) so that their values are of type
    'float' (if possible).
    """
    model = _format_notes_wrapper(model, "model")
    model = _format_notes_wrapper(model, "metabolites")
    model = _format_notes_wrapper(model, "reactions")
    return model


def _format_notes_wrapper(model, notes_on="model"):
    """
    Formats the dictionary stored in 'notes' of a given Cobra objects
    so that their values are of type 'float' (if possible).
    """
    if notes_on == "model":
        to_be_populated = [model]
    else:
        try:
            to_be_populated = getattr(model, notes_on)
        except AttributeError:
            to_be_populated = None

    allowable_param_keys = _generate_allowable_param_keys(model)

    if _check_is_iterable(to_be_populated):
        for x in to_be_populated:
            if _check_has_notes(x):

                populate_attributes(x, allowable_param_keys)

    return model


def _check_is_iterable(obj):
    """
    Returns True if the object is iterable; returns False otherwise.
    """
    try:
        getattr(obj, "__iter__")
        out = True
    except AttributeError as e:
        print(f"{e}: object {type(obj)} is not iterable.")
        out = False
    return out


def _check_has_notes(obj):
    """
    Returns True if the object contains an attribute 'notes'; returns
    False otherwise.
    """
    try:
        getattr(obj, "notes")
        out = True
    except AttributeError as e:
        print(f"{e}: object {type(obj)} does not have attribute 'notes'.")
        out = False
    return out


def populate_attributes(obj, allowable_param_keys):
    for parameter, value in obj.notes.items():

        # Get the appropriate value to be stored:
        try:
            magnitude, units = value.split(" ", 1)
            value_to_store = Q_(float(magnitude), units)
        except ValueError:
            value_to_store = value
        except AttributeError:
            value_to_store = float(value)

        # Check if scalar or indexed parameter
        key = None
        suffixes = tuple(allowable_param_keys)
        if parameter.endswith(suffixes):
            parameter, key = parameter.split("_", 1)

        # Store the values
        try:
            getattr(obj, parameter)
            _helper_populate_attr(obj, "_"+parameter, value_to_store, key)
        except AttributeError as e:
            print(e)

    return


def _helper_populate_attr(obj, attr_name, value_to_store, key):
    if key is None:
        to_store = value_to_store
    else:
        to_store = {key: value_to_store}

    if isinstance(getattr(obj, attr_name), dict):
        getattr(obj, attr_name).update(to_store)
    else:
        setattr(obj, attr_name, to_store)


def _generate_allowable_param_keys(model):
    compartments = [x for x in model.compartments]
    compartment_combinations = ["".join(x) for x in itertools.permutations(compartments, 2)]
    suffixes = compartments + compartment_combinations
    return ["_" + x for x in suffixes]


# LEGACY function
def _convert_entries_to_float(dict_in):
    """
    Converts a dictionary's values into floats, wherever possible.
    """
    dict_out = {}
    for k, val in dict_in.items():
        try:
            dict_out[k] = float(val)
        except ValueError:
            dict_out[k] = val

    return dict_out