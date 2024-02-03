import itsinv
import xarray
import numpy as np

# Invert for velocity
xy = (-3310000, 257700)
t0 = "2015-10-01"
t1 = "2023-01-01"

solu, res, xsub, mask = itsinv.grid_inversion_ncpu(
    xy,
    half_dist=6000,
    lt=1e-4,
    lx=10,
    sat_filt=["1A", "1B", "4", "5", "7", "8", "9"],
    start_date=t0,
    stop_date=t1,
    pbar=True,
    return_data=True,
    ncpu=4
)

# Save itslive datacube
xsub.to_netcdf("hubbard_itslive_%s_%s.nc" % (t0, t1))

# Convert decimal year to datetime64
s_in_year = (60*60*24*365)
solu_datetime = [None]*len(solu)
for i in range(len(solu)):
    whole_year = int(solu[i])
    frac_year = solu[i] - whole_year
    solu_datetime[i] = np.datetime64(str(whole_year))
    solu_datetime[i] += np.timedelta64(int(frac_year*s_in_year), "s")

xres = xarray.Dataset(
    data_vars = dict(
        vx=(["time", "y", "x"], res["vx"]),
        vy=(["time", "y", "x"], res["vy"]),
    ),
    coords = dict(
        x = xsub.x,
        y = xsub.y,
        time = solu_datetime,
    ),
    attrs = dict(
        projection = xsub.projection,
        GDAL_AREA_OR_POINT = xsub.GDAL_AREA_OR_POINT,
    ))

xres.to_netcdf("hubbard_inversion_%s_%s.nc" % (t0, t1))
