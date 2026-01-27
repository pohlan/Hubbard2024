import numpy as np
import datetime as dt
import pandas as pd
import geopandas as gpd
import rasterio
import matplotlib.pyplot as plt
from rasterio.plot import show


terminus_data_path = "/home/annegret/Projects/Hubbard2024/data/terminus_data/"        # folder where all the data for terminus position and terminus sections is stored
veloc_data_path    = "/home/annegret/Projects/Hubbard2024/data/velocity/"
image_path         = "/home/annegret/Projects/Hubbard2024/data/images/"
fig_path           = "/home/annegret/Projects/Hubbard2024/terminus_change/figures/"   # folder to save figures

# load sections, terminus positions and velocity data
df_veloc      = pd.read_csv(veloc_data_path+"velocity_by_section.csv")
df_veloc_ortho = pd.read_csv(veloc_data_path+"velocity_ortho_by_section.csv")
df_terminus   = pd.read_csv(terminus_data_path+"advance_by_section.csv")
df_ablation   = pd.read_csv(terminus_data_path+"frontal_ablation.csv")
df_retr_rate  = pd.read_csv(terminus_data_path+"advance_rate_by_section.csv")
df_veloc.Date = pd.to_datetime(df_veloc.Date)

# gpkg of sections to plot map
gdf_sections = gpd.read_file(terminus_data_path+"terminus_sections_rectangle.gpkg")
src          = rasterio.open(image_path+"composite.tif")

# loop through sections, plot and save velocity time series
fig1, ax1 = plt.subplots(5,1, layout='constrained', figsize=(9,9))
fig2, ax2 = plt.subplots(5,1, layout='constrained', figsize=(9,9))
fig3, ax3 = plt.subplots(5,1, layout='constrained', figsize=(9,9))

i = 0
for (sterm, sveloc, sretr, sabl) in zip(df_terminus, df_veloc, df_retr_rate, df_ablation):
    if sterm == "Date" or sveloc == "Date":
        continue

    # plot velocity with advance
    cols = ["royalblue", "black"]
    ax1[i].plot(df_veloc.Date, df_veloc[sveloc], color=cols[0], label="Median velocity v")
    ax1[i].set_ylabel(r"v (m a$^{-1}$)", color=cols[0])
    ax1[i].set_title(sveloc.split(" ")[0].replace("_", " "))
    secax = ax1[i].twinx()
    secax.plot(df_veloc.Date, df_terminus[sterm], color=cols[1], label="Advance position (m)")
    secax.set_ylabel("Terminus advance (m)", color=cols[1])
    secax.set_xlim(np.min(df_veloc.Date), np.max(df_veloc.Date))

    # plot velocity with time derivative of advance time series
    ax2[i].plot(df_veloc.Date, df_veloc[sveloc], color=cols[0], label="Median velocity v")
    ax2[i].plot(df_veloc.Date[1:-1], df_retr_rate[sretr], color=cols[1], label="Advance rate")
    ax2[i].set_ylabel(r"m a$^{-1}$")
    ax2[i].set_title(sveloc.split(" ")[0].replace("_", " "))
    ax2[i].legend(loc="lower left")

    ax3[i].plot(df_veloc.Date, df_veloc[sveloc], color=cols[0], label="Median velocity v")
    ax3[i].plot(df_veloc.Date[1:-1], df_ablation[sabl], color=cols[1], label="Frontal ablation")
    ax3[i].set_ylabel(r"m a$^{-1}$")
    ax3[i].set_title(sveloc.split(" ")[0].replace("_", " "))
    ax3[i].legend(loc="lower left")

    i+=1

fig1.savefig(fig_path+"veloc_vs_terminus_advance.pdf")
fig2.savefig(fig_path+"veloc_vs_terminus_advance_rate.pdf")
fig3.savefig(fig_path+"veloc_vs_frontal_ablation.pdf")



# plot map of sections
fig, ax = plt.subplots(figsize=(10, 10))

show(src, ax=ax)
gdf_sections.plot(ax=ax, color='skyblue', edgecolor='black', alpha=0.6)

dx = [-600,-500,-200,100,0]
dy = [0,-500,-500,-500,-100]

for idx, row in gdf_sections.iterrows():
    # Get the center point of the polygon
    centroid = row.geometry.centroid

    # Annotate using a column name (e.g., 'name' or 'title')
    ax.annotate(text=row['section'].replace("_"," "),
                xy=(centroid.x+dx[idx], centroid.y+dy[idx]),
                ha='center', # fontweight='bold',
                fontsize=15)

# add a scale bar (without solid_capstyle the lines overlap in x-direction)
ax.plot([580000, 581000], [6650020, 6650020], 'w', solid_capstyle='butt', linewidth=4)
ax.plot([581000, 582000], [6650020, 6650020], 'k', solid_capstyle='butt', linewidth=4)
ax.plot([582000, 584000], [6650020, 6650020], 'w', solid_capstyle='butt', linewidth=4)
ax.plot([580000, 581000], [6650100, 6650100], 'k', solid_capstyle='butt', linewidth=4)
ax.plot([581000, 582000], [6650100, 6650100], 'w', solid_capstyle='butt', linewidth=4)
ax.plot([582000, 584000], [6650100, 6650100], 'k', solid_capstyle='butt', linewidth=4)

ax.text(580000, 6650300, '0', fontsize=14, color='w', fontweight='bold', ha='center')
ax.text(581000, 6650300, '1', fontsize=14, color='w', fontweight='bold', ha='center')
ax.text(582000, 6650300, '2', fontsize=14, color='w', fontweight='bold', ha='center')
ax.text(584000, 6650300, '4', fontsize=14, color='w', fontweight='bold', ha='center')
ax.text(584100, 6650020, 'km', fontsize=14, color='w', fontweight='bold', va='center')

ax.axis('off')

ax.set_xlim(5.79e5, 5.9e5)
ax.set_ylim(6.649e6, 6.659e6)

plt.savefig(fig_path+"section_map.pdf")
