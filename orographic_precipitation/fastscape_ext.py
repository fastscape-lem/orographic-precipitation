import numpy as np
import xsimlab as xs

from fastscape.models import basic_model
from fastscape.processes import FlowAccumulator, SurfaceTopography, UniformRectilinearGrid2D
from .orographic_precipitation import compute_orographic_precip


@xs.process
class OrographicPrecipitation:
    """Computes orographic precipitation following Smith & Barstad (2004)
    """
    # --- initial conditions
    lapse_rate = xs.variable(description="adiabatic lapse rate")
    lapse_rate_m = xs.variable(description="moist adiabatic lapse rate")
    ref_density = xs.variable(description="reference density")

    # --- input variables
    lat = xs.variable(description="latitude")
    precip_base = xs.variable(description="background precipitation")
    wind_speed = xs.variable(description="windspeed")
    wind_dir = xs.variable(description="wind direction")
    conv_time = xs.variable(description="conversion time")
    fall_time = xs.variable(description="fallout time")
    nm = xs.variable(description="moist stability frequency")
    hw = xs.variable(description="water vapor scale height")
    cw = xs.variable(description="uplift sensitivity", intent="out")

    # --- variables needed for computation
    dx = xs.foreign(UniformRectilinearGrid2D, "dx")
    dy = xs.foreign(UniformRectilinearGrid2D, "dy")
    shape = xs.foreign(UniformRectilinearGrid2D, "shape")
    elevation = xs.foreign(SurfaceTopography, "elevation")

    # --- output variable
    precip_rate = xs.variable(
        dims=("y", "x"), description="precipitation rate", intent="out", attrs={"units": "mm/h"}
    )

    def _get_params(self):
        return {
            "latitude" : self.lat,
            "p0" : self.precip_base,
            "windspeed" : self.wind_speed,
            "winddir" : self.wind_dir,
            "tau_c" : self.conv_time,
            "tau_f" : self.fall_time,
            "nm" : self.nm,
            "hw" : self.hw,
            "cw" : self.ref_density * self.lapse_rate_m / self.lapse_rate}

    def initialize(self):
        self._params = self._get_params()
        self.precip_rate = np.zeros(self.shape)

    def run_step(self):
        self._params.update(self._get_params())
        self.precip_rate = compute_orographic_precip(self.elevation,
                                                self.dx,
                                                self.dy,
                                                **self._params)


@xs.process
class OrographicDrainageDischarge(FlowAccumulator):
    """Flowaccumulation including orographic precipitation
    """
    runoff = xs.foreign(OrographicPrecipitation, 'precip_rate')
    discharge = xs.on_demand(
        dims=('y','x'), description='discharge from orographic precipitation', attrs={"units": "mm/h"}
    )

    def run_step(self):
        super().run_step()

        # scale mm/h to m/yr
        self.flowacc *= 8.76

    @discharge.compute
    def _discharge(self):
        return self.flowacc / 8.76


precip_model = basic_model.update_processes({
    'orographic': OrographicPrecipitation,
    'drainage': OrographicDrainageDischarge
})