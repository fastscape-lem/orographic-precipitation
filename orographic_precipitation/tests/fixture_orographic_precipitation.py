import pytest

@pytest.fixture
def input_params():
    lapse_rate = -5.8e-3
    lapse_rate_m = -6.5e-3
    ref_density = 7.4e-3

    params = {
        'latitude': 0,
        'precip_base': 0,
        'wind_speed': 15,
        'wind_dir': 270,
        'conv_time': 1000,
        'fall_time': 1000,
        'nm': 0.005,
        'hw': 2500,
        'lapse_rate': lapse_rate,
        'lapse_rate_m': lapse_rate_m,
        'ref_density': ref_density,
        'cw': ref_density * lapse_rate_m / lapse_rate
        }

    return params
