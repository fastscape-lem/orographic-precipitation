import pytest

@pytest.fixture
def input_params():
    gamma = -5.8e-3    #-6.49
    Gamma_m = -6.5e-3  #-5
    rhosref = 7.4e-3

    params = {
        'latitude': 0,
        'p0': 0,
        'windspeed': 15,
        'winddir': 270,
        'tau_c': 1000,
        'tau_f': 1000,
        'nm': 0.005,
        'hw': 2500,
        'cw': rhosref * Gamma_m / gamma
        }

    return params
