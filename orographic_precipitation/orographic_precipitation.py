
#!/usr/bin/env python

import numpy as np

EPS = np.finfo(float).eps

def compute_orographic_precip(elevation, dx, dy, **param):
    """Compute orographic precipitation.

    Parameters
    ----------
    elevation : array_like
        2D input array of a given elevation
    dx, dy : int
        Horizontal and vertical resolution in [m]
    **param
        A dictionary used to store relevant parameters for computation.

    param kwargs
    ----------------
    latitude (float) : Coriolis effect decreases as latitude decreases
    p0 (float) : uniform precipitation rate [mm hr-1], usually [0, 10]
    windspeed (float) : [m s-1]
    winddir (float) : wind direction [0: north, 270: west]
    tau_c (float) : conversion time delay [s]
    tau_f (float) : fallout time delay [s]
    nm (float) : moist stability frequency [s-1]
    hw (float) : water vapor scale height [m]
    cw (float) : uplift sensitivity [kg m-3], product of saturation water vapor sensitivity rhosref [kg m-3] and environmental lapse rate (gamma/gamma_n)

    Returns
    -------
    array_like
        2D array structure the same size as elevation with precipitation [mm hr-1]
    """

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
