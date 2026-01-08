import numpy as np
import xarray

import itsinv

# Invert for velocity
xy = (-3310000, 257700)
t0 = "2015-10-01"
t1 = "2023-01-01"

xinv, xsub, mask = itsinv.grid_inversion(
    xy,
    half_dist=8000,
    lt=1e-4,
    lx=0,
    sat_filt=["1A"],
    start_date=t0,
    stop_date=t1,
    pbar=True,
    return_data=True,
    ncpu=4,
)

# Save itslive datacube
xsub.to_netcdf("hubbard_itslive_%s_%s.nc" % (t0, t1))

# Save inverted datacube11
xinv.to_netcdf("hubbard_inversion_%s_%s.nc" % (t0, t1))
