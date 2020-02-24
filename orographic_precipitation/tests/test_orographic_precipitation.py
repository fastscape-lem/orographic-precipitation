import numpy as np

import pytest

from orographic_precipitation import compute_orographic_precip
from fixture_orographic_precipitation import input_params

def test_compute_orographic_precip():
    input = input_params()

    assert type(input).__module__ == np.__name__
