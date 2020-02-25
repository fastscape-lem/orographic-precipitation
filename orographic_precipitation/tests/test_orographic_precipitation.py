import numpy as np

import pytest

from orographic_precipitation import compute_orographic_precip
from orographic_precipitation.tests.fixture_orographic_precipitation import input_params

#pytest.mark.usefixtures(input_params)
def test_compute_orographic_precip(input_params):
    dx = dy = 0.01
    length = 1.5
    x = np.arange(0.1, length, dx)
    y = np.arange(0.1, length, dy)
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)

    output = compute_orographic_precip(z, dx, dy, **input_params)

    assert type(output).__module__ == np.__name__
    assert z.shape == output.shape

if __name__ == '__main__':
    main()
