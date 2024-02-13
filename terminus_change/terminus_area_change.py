import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.ops import split
from shapely import Point
import datetime as dt
import statsmodels.api as sm

local_path = "/home/annegret/myfolder/Hubbard2024/terminus_change/"

# read in terminus traces (UTM07)
terminus_path = local_path + "terminus_data/Hubbard_2017_2021_UTM07.shp"
df_terminus = gpd.read_file(terminus_path)   # df_terminus.geometry is already in LineString format
# Mcnabb terminus data
terminus_mcnabb_path = local_path + "terminus_data/hubbard.shp"
df_terminus_mcnabb = gpd.read_file(terminus_mcnabb_path)

# read in section polygons
path_sections = local_path + "terminus_data/terminus_sections_rectangle.gpkg"
df_sec = gpd.read_file(path_sections) # df_sA.geometry is already in Polygon format
df_sec = df_sec.to_crs(df_terminus.crs)

def area_rule(name, geom1, geom2):
    x1, y1 = geom1.exterior.coords.xy
    x2, y2 = geom2.exterior.coords.xy

    if name=="section_A" or name=="section_B1" or name=="section_B2":
        xmax1 = np.max(x1)
        xmax2 = np.max(x2)
        if xmax1 > xmax2:
            return geom1
        else:
            return geom2
    else:
        ymax1 = np.max(y1)
        ymax2 = np.max(y2)
        if ymax1 > ymax2:
            return geom1
        else:
            return geom2

def get_width(geom):
    poly = list(geom.geoms)  # make Polygons out of MultiPolygon to be able to extract coordinates
    x, y = poly[0].exterior.coords.xy
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
    width = min(edge_length)
    return width

# split polygon with terminus line
list_areas = []
n = 0
for (day, terminus_line) in zip(df_terminus.Date, df_terminus.geometry):
    newrow = {}
    newrow["Date"] = dt.datetime.strptime(day, "%Y-%m-%d").date()
    for (sect_name, sect_geom) in zip(df_sec.section, df_sec.geometry):
        split_pols = split(sect_geom, terminus_line)
        if len(list(split_pols.geoms)) != 2:
            continue
        # assert list(split_pols.geoms)[0].area > list(split_pols.geoms)[1].area
        glacier_pol = area_rule(sect_name, split_pols.geoms[0], split_pols.geoms[1])
        width       = get_width(sect_geom)
        newrow[sect_name] = glacier_pol.area / width
    list_areas.append(newrow)
    n += 1
df_areas = pd.DataFrame(list_areas)
df_areas = df_areas.sort_values(by=["Date"])

# plt.plot(df_areas.Date,df_areas.section_A / np.max(df_areas.section_A), label="section A")
# plt.plot(df_areas.Date,df_areas.sectionB_2 / np.max(df_areas.sectionB_2), label="section B2")
# plt.plot(df_areas.Date,df_areas.sectionD / np.max(df_areas.sectionD), label="section D")
# plt.legend()
# plt.savefig("terminus_vs_time.jpg")

# use lowess filter to smooth time series
lowess = sm.nonparametric.lowess
t0 = dt.datetime.strptime("2017-01-01", "%Y-%m-%d").date()
tend = dt.datetime.strptime("2021-12-31", "%Y-%m-%d").date()

t_index = np.zeros(len(df_areas.Date))
for (i, da) in enumerate(df_areas.Date):
    deltat = da - t0
    t_index[i] = deltat.days
t_index

fig, ax = plt.subplots(2, 1, sharey=True)
for s in df_sec.section:
    z = lowess(df_areas[s] - np.mean(df_areas[s]), t_index, frac=1/20)
    t_smooth = np.zeros(z.shape[0], dtype=dt.date)
    for (i,zi) in enumerate(z[:,0]):
        t_smooth[i] = t0 + dt.timedelta(days=zi)
    # plt.figure()
    # plt.plot(df_areas.Date,df_areas[s] - np.mean(df_areas[s]), label=s, lw=1, alpha=0.9)
    if s == "section_A" or s == "section_B1":
        ax[0].plot(t_smooth, z[:,1], label=s)
        ax[0].legend()
    elif  s == "section_B2" or s == "section_C":
        ax[1].plot(t_smooth, z[:,1], label=s)
        ax[1].legend()

ax[0].set_ylabel("rel. terminus pos. [m]")
ax[1].set_ylabel("rel. terminus pos. [m]")
ax[1].set_xlabel("Date")
plt.savefig(local_path + "smoothed_termini.jpg")
