# Add geotransform to inverted velocity file
import rioxarray
import xarray
import affine
import numpy as np

xds = xarray.open_dataset(
    "hubbard_inversion_2015-10-01_2023-01-01.nc", decode_coords="all"
)

crs = "epsg:3414"
xform = affine.Affine.from_gdal(
    np.min(xds.x.to_numpy()),
    np.diff(xds.x)[0],
    0,
    np.max(xds.y.to_numpy()),
    0,
    np.diff(xds.y)[0],
)

xds = xds.rio.write_crs("epsg:3413")
xds = xds.rio.write_transform(xform)

xds.to_netcdf("hubbard_inversion_2015-10-01_2023-01-01_xform.nc")
