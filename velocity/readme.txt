Michael Christoffersen

To make a movie:
First run make_velocity_datacube to make an inverted datacube of velocities.
Then run make_movie_frames to generate all of the frames of the movie. The extent of the plotting window and the projection are hard-coded in the script so you'll have to go in and modify it if not plotting Hubbard.
Then run the ffmpeg command in ffmpeg.bash, which will join all of the frames into a movie. You can control the speed of the movie with the frames per second parameter in ffmpeg.


make_velocity_datacube.py - Downloads a square of ITS_LIVE data and inverts for a velocity timeseries at every pixel in the square. Outputs a netcdf xarray.

make_movie_frames.py - Takes the datacube output by make_velocity_datacube and generates daily velocity images.

ffmpeg.bash - Example ffmpeg command to make a movie from a bunch of images.

invert_annual_velocity.ipynb - Very similar to make_velocity_datacube but in a notebook.

double_peak_plot.ipynb - Make plots of amplitude and timing of double velocity peaks.

itsinv.py - Velocity inversion routines. Includes functions to download itslive data and invert for velocity at a point or on a grid.

datacube_tools.py is necessary for the inversion code (it downloads itslive data) and is from here:
https://github.com/nasa-jpl/its_live/blob/main/notebooks/datacube_tools.py

pod.ipynb - Something like the POD paper. Except not using SVD to interpolate to constant times and not smoothing V rows

add_transform.py - Add geotransform info to datacube
