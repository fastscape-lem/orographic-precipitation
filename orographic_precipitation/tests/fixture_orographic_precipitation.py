import pytest
import xarray as xr
import numpy as np


@pytest.fixture
def input_params():
    lapse_rate = -5.8e-3
    lapse_rate_m = -6.5e-3
    ref_density = 7.4e-3

    params = {
        'latitude': 0,
        'precip_base': np.ones((140,140)),
        'precip_min': 0.01,
        'rainfall_frequency': 0.1,
        'wind_speed': 15,
        'wind_dir': 270,
        'conv_time': 1000,
        'fall_time': 1000,
        'nm': 0.005,
        'hw': 2500,
        'lapse_rate': lapse_rate,
        'lapse_rate_m': lapse_rate_m,
        'ref_density': ref_density,
        'cw': ref_density * lapse_rate_m / lapse_rate,
        }

    return params
