import numpy as np
import xsimlab as xs

from fastscape.processes import SurfaceTopography, UniformRectilinearGrid2D
from orographic_precipitation import compute_orographic_precip


@xs.process
class OrographicPrecipitation:
    """Computes orographic precipitation following Smith & Barstad (2004)
    """
    # --- initial conditions
    gamma = xs.variable(description="adiabatic lapse rate")
    gamma_m = xs.variable(description="moist adiabatic lapse rate")
    rhosref = xs.variable(description="reference density")
    
    # --- input variables
    latitude = xs.variable(description="latitude")
    p0 = xs.variable(description="background precipitation")
    windspeed = xs.variable(description="windspeed")
    winddir = xs.variable(description="wind direction")
    tau_c = xs.variable(description="conversion time")
    tau_f = xs.variable(description="fallout time")
    nm = xs.variable(description="moist stability frequency")
    hw = xs.variable(description="water vapor scale height")
    cw = xs.variable(description="uplift sensitivity", intent="out")

    # --- variables needed for computation
    dx = xs.foreign(UniformRectilinearGrid2D, "dx")
    dy = xs.foreign(UniformRectilinearGrid2D, "dy")
    shape = xs.foreign(UniformRectilinearGrid2D, "shape")
    elevation = xs.foreign(SurfaceTopography, "elevation")
    
    # --- output variable
    precip_rate = xs.variable(dims=("y", "x"), description="precipitation rate", intent="out")
    precip = xs.variable(dims=("y", "x"), description="precipitation", intent="out")
    
    def initialize(self):
        self.cw = self.rhosref * self.gamma_m / self.gamma
        self.precip_rate = np.zeros(self.shape)
    
    @xs.runtime(args="step_delta")
    def run_step(self, dt):
        _params = {"latitude" : self.latitude, 
                   "p0" : self.p0, 
                   "windspeed" : self.windspeed, 
                   "winddir" : self.winddir, 
                   "tau_c" : self.tau_c, 
                   "tau_f" : self.tau_f, 
                   "nm" : self.nm, 
                   "hw" : self.hw, 
                   "cw" : self.cw}
        self.precip_rate = compute_orographic_precip(self.elevation,
                                                self.dx,
                                                self.dy,
                                                **_params)
        self.precip = self.precip_rate[:] * dt