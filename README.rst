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

Usage
-----

1. Load relevant packages and initialize functions to compute orographic
   precipitation, set up a topography and plot the resulting precipitation matrix.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    EPS = np.finfo(float).eps


    def compute_orographic_precip(elevation, dx, dy, **param):
        """Compute orographic precipitation."""
        # --- wind components
        u0 = -np.sin(param['winddir'] * 2 * np.pi / 360) * param['windspeed']
        v0 = np.cos(param['winddir'] * 2 * np.pi / 360) * param['windspeed']

        # --- other factors
        f_coriolis = 2 * 7.2921e-5 * np.sin(param['latitude'] * np.pi / 180)

        # --- pad raster boundaries prior to FFT
        calc_pad = int(np.ceil(((sum(elevation.shape))) / 2) / 100 * 100)
        pad = min([calc_pad, 200])

        h = np.pad(elevation, pad, 'constant')
        nx, ny = h.shape

        # --- FFT
        hhat = np.fft.fft2(h)

        x_n_value = np.fft.fftfreq(ny, (1. / ny))
        y_n_value = np.fft.fftfreq(nx, (1. / nx))

        x_len = nx * dx
        y_len = ny * dy
        kx_line = 2 * np.pi * x_n_value / x_len
        ky_line = 2 * np.pi * y_n_value / y_len
        kx = np.tile(kx_line, (nx, 1))
        ky = np.tile(ky_line[:, None], (1, ny))

        # --- vertical wave number (m)
        sigma = kx * u0 + ky * v0

        mf_num = param['nm']**2 - sigma**2
        mf_den = sigma**2 - f_coriolis**2

        # numerical stability
        mf_num[mf_num < 0] = 0.
        mf_den[(mf_den < EPS) & (mf_den >= 0)] = EPS
        mf_den[(mf_den > -EPS) & (mf_den < 0)] = -EPS
        sign = np.where(sigma >= 0, 1, -1)

        m = sign * np.sqrt(np.abs(mf_num / mf_den * (kx**2 + ky**2)))

        # --- transfer function
        P_karot = ((param['cw'] * 1j * sigma * hhat) /
                   ((1 - (param['hw'] * m * 1j)) *
                    (1 + (sigma * param['tau_c'] * 1j)) *
                    (1 + (sigma * param['tau_f'] * 1j))))

        # --- inverse FFT, de-pad, convert units, add uniform rate
        P = np.fft.ifft2(P_karot)
        P = np.real(P[pad:-pad, pad:-pad])
        P *= 3600   # mm hr-1
        P += param['p0']
        P[P < 0] = 0

        return P

    def gauss_topography(dx, dy):
        h_max = 500.
        x0 = -25e3
        y0 = 0
        sigma_x = sigma_y = 15e3

        x, y = np.arange(-100e3, 200e3, dx), np.arange(-150e3, 150e3, dy)
        X, Y = np.meshgrid(x, y)

        h = (h_max *
             np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                      ((Y - y0)**2 / (2 * sigma_y**2)))))

        return X, Y, h

    def plot2d(X, Y, field):
        fig, ax = plt.subplots(figsize=(6, 6))
        pc = ax.pcolormesh(X, Y, field)
        ax.set_aspect(1)
        fig.colorbar(pc)

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
