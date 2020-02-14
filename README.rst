Linear Theory of Orographic Precipitation
=========================================

A `Python`_ framework that implements the Linear Theory of Orographic Precipitation
following `Smith & Barstad (2004)`_.

.. _`Python`: https://www.python.org
.. _`Smith & Barstad (2004)`: https://journals.ametsoc.org/doi/full/10.1175/1520-0469%282004%29061%3C1377%3AALTOOP%3E2.0.CO%3B2

The model includes airflow dynamics, condensed water advection, and downslope
evaporation. It consists of two vertically-integrated steady-state advection
equations describing: (i) the cloud water density and (ii) the hydrometeor
density. Solving these equations using Fourier transform techniques,
derives a single formula relating terrain and precipitation.

Please refer to the original manuscript of Smith and Barstad (2004) to understand
the model physics and limitations.

Installation
------------

Required dependencies:
* Python 3.6 or later
* `numpy`_
* `matplotlib`_

.. _`numpy`: https://numpy.org
.. _`matplotlib`: https://matplotlib.org/3.1.1/index.html

orographic precipitation can be installed using `pip`

.. code-block::

  $ pip install orographic_precipitation

Usage
-----

1. Import relevant functions to compute orographic precipitation,
   set up a topography and plot the resulting precipitation matrix.

.. code-block:: python

    from orographic_precipitation import (compute_orographic_precip,
                                         gauss_topography,
                                         plot2d)

2. Create example topography, *i.e.* an isolated circular Gaussian hill.

.. code-block:: python

    dx = 750.
    dy = 750.

    X, Y, elevation = gauss_topography(dx, dy)

    plot2d(X, Y, elevation)

.. image:: doc/_static/gauss_topo.png
   :width: 400px

3. Initialize dictionary with relevant parameters, compute and plot orographic
   precipitation.

.. code-block:: python

    gamma = -5.8    #-6.49
    Gamma_m = -6.5  #-5
    rhosref = 7.4e-3

    param = {
    'latitude': 40,
    'p0': 7,                          # uniform precipitation rate
    'windspeed': 10,
    'winddir': 270,
    'tau_c': 1000,                    # conversion time
    'tau_f': 1000,                    # fallout time
    'nm': 0.005,                      # moist_stability_freq
    'hw': 5000,                       # water_vapor_scale_height
    'cw': rhosref * Gamma_m / gamma   # uplift_sensitivity
    }

    P = compute_orographic_precip(elevation, dx, dy, **param)

    plot2d(X, Y, P)

.. image:: doc/_static/orograph_precip_example.png
   :width: 400px

Acknowledgement
---------------

This project is supported by the `Earth Surface Process Modelling`_ group at
the German Geoscience Research Institute in Potsdam, Germany.

.. _`Earth Surface Process Modelling`: http://www.gfz-potsdam.de/en/section/earth-surface-process-modelling/
