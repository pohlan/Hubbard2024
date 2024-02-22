import geopandas as gpd
import numpy as np
import datetime as dt
import pandas as pd
import rioxarray
import matplotlib.pyplot as plt

terminus_data_path = "/home/annegret/Projects/Hubbard2024/data/terminus_data/"        # folder where all the data for terminus position and terminus sections is stored
veloc_data_path    = "/home/annegret/Projects/Hubbard2024/data/velocity/"
fig_path           = "/home/annegret/Projects/Hubbard2024/terminus_change/figures/"   # folder to save figures

# load sections, terminus positions and velocity data
geodf_sections = gpd.read_file(terminus_data_path+"section_polygons_glacier_only.gpkg")
xds_veloc      = rioxarray.open_rasterio(veloc_data_path+"hubbard_inversion_2015-10-01_2023-01-01_xform.nc", masked=True)
df_terminus    = pd.read_csv("data/terminus_data/terminus_position_smooth.csv")

# loop through sections, plot and save velocity time series
fig, ax = plt.subplots(5,1, layout='constrained', figsize=(9,9))
t_term = [dt.datetime.strptime(d, "%Y-%m-%d").date() for d in df_terminus["Date"]]                     # time series of terminus positions
v_t    = [dt.datetime.strptime(t.isoformat(), "%Y-%m-%dT%H:%M:%S") for t in xds_veloc["time"].values]  # time series of velocities
dict_vel = {"Date":v_t}
for (i,(sec,geom)) in enumerate(zip(geodf_sections.section, geodf_sections.geometry.values)):
    clipped = xds_veloc.rio.clip([geom], geodf_sections.crs)     # extract velocity data points that lie within the section
    v = np.sqrt(clipped.vx**2 + clipped.vy**2)                   # get absolute velocity
    v_median = v.median(dim=["x","y"])
    dict_vel[sec+" median |v| [m/yr]"] = v_median.values         # save median velocity values in dictionary to save later

    # plot
    cols = ["darkblue", "orange"]
    ax[i].plot(v_t, v_median, color=cols[0])
    ax[i].set_ylabel("median |v| [m/yr]", color=cols[0])
    ax[i].set_title(sec)
    secax = ax[i].twinx()
    secax.plot(t_term, df_terminus[sec+" [m]"], color=cols[1])
    secax.set_ylabel("rel. terminus pos. [m]", color=cols[1])
    secax.set_xlim(np.min(t_term), np.max(t_term))
plt.savefig(fig_path+"veloc_vs_terminus.jpg")
