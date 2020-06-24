import numpy as np
import pytest

from orographic_precipitation import compute_orographic_precip
from orographic_precipitation.fastscape_ext import OrographicPrecipitation
from orographic_precipitation.tests.fixture_orographic_precipitation import input_params

def test_orographic_precipitation(input_params):

    process_params = input_params.copy()
    process_params.pop("cw")
    grid = (3, 2)
    elevation = np.random.uniform(size=grid)
    dx = dy = 0.1
    
    p = OrographicPrecipitation(
        shape = grid,
        dx = dx,
        dy = dy,
        elevation = elevation,
        **process_params
    )
    
    p.initialize()

    np.testing.assert_equal(p.precip_rate.shape, grid)

    p.run_step()

    expected = compute_orographic_precip(elevation, dx, dy, **input_params)
    np.testing.assert_equal(p.precip_rate, expected)
    assert p._get_params().items() <= input_params.items()
