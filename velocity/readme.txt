Michael Christoffersen

make_velocity_datacube.py - Downloads a square of ITS_LIVE data and inverts for a velocity timeseries at every pixel in the square. Outputs a netcdf xarray.

invert_annual_velocity.ipynb - Very similar to make_velocity_datacube but in a notebook.

double_peak_plot.ipynb - Make plots of amplitude and timing of double velocity peaks.

itsinv.py - Velocity inversion routines. Includes functions to download itslive data and invert for velocity at a point or on a grid.

datacube_tools.py is necessary for the inversion code (it downloads itslive data) and is from here:
https://github.com/nasa-jpl/its_live/blob/main/notebooks/datacube_tools.py
