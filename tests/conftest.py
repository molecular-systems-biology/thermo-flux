import pytest
import os
from cobra.io import read_sbml_model
from thermo_flux.core.model import ThermoModel
from equilibrator_api import Q_

@pytest.fixture(scope="session")
def tmodel():
    """Load the default tmodel once per test session."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(test_dir, 'yeast_merged.xml')
    
    model = read_sbml_model(model_path)
    return ThermoModel(model, 
                      pH={'c': Q_(7), 'm': Q_(7.4), 'e': Q_(5)},
                      I={'c': Q_(0.25, 'M'), 'm': Q_(0.25, 'M'), 
                         'e': Q_(0.25, 'M')}, 
                      T=Q_(303.15, 'K'), 
                      pMg={'c': Q_(3), 'm': Q_(3), 'e': Q_(3)},
                      phi={'ec': Q_(-0.06, 'V'),
                           'cm': Q_(-0.16, 'V')}, 
                      update_thermo_info=True)

@pytest.fixture(scope="session")
def tmodel_b():
    """Load a different model once per test session."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(test_dir, 'yeast_merged.xml')
    
    model = read_sbml_model(model_path)
    return ThermoModel(model,
                      pH={'c': Q_(7), 'm': Q_(1), 'e': Q_(5)},
                      I={'c': Q_(0.25, 'M'), 'm': Q_(0.25, 'M'), 'e': Q_(0.25, 'M')},
                      T=Q_(303.15, 'K'),
                      pMg={'c': Q_(3), 'm': Q_(3), 'e': Q_(3)},
                      phi={'ec': Q_(-0.06, 'V'), 'cm': Q_(-0.16, 'V')},
                      update_thermo_info=False)

@pytest.fixture(scope="session")
def no_model():
    """Fixture that does not load any model."""
    return None