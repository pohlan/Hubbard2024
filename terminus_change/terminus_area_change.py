import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.ops import split
from shapely import Point
import datetime as dt
import statsmodels.api as sm

data_path = "/home/annegret/Projects/Hubbard2024/data/terminus_data/"     # folder where all the data for terminus position and terminus sections is stored
fig_path  = "/home/annegret/Projects/Hubbard2024/terminus_change/figures/"           # folder to save figures

# read in terminus traces (UTM07)
terminus_path = data_path + "Hubbard_2017_2021_UTM07.gpkg"
df_terminus = gpd.read_file(terminus_path)   # df_terminus.geometry is already in LineString format
# Mcnabb terminus data
terminus_mcnabb_path = data_path + "Hubbard_McNabb.gpkg"
df_terminus_mcnabb = gpd.read_file(terminus_mcnabb_path)

# read in section polygons
path_sections = data_path + "terminus_sections_rectangle.gpkg"
df_sec = gpd.read_file(path_sections) # df_sA.geometry is already in Polygon format
df_sec = df_sec.to_crs(df_terminus.crs)
df_sec.sort_values(by=["section"], inplace=True)

# define some functions for deriving the area change
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

def append_positions_by_section(newrow, terminus_line, list_small_polygons):
    for (sect_name, sect_geom) in zip(df_sec.section, df_sec.geometry):
        split_pols = split(sect_geom, terminus_line)
        if len(list(split_pols.geoms)) != 2:
            continue
        glacier_pol = area_rule(sect_name, split_pols.geoms[0], split_pols.geoms[1])
        width       = get_width(sect_geom)
        newrow[sect_name] = glacier_pol.area / width
        # for one terminus position at time step 102, which is one of the furthest retreated ones:
        # save the glacier part of the section as a shape file
        sctions = [l.get("section") for l in list_small_polygons]
        if day == df_terminus.Date[102] and not sect_name in sctions :
            polrow = {"section":sect_name, "geometry": glacier_pol}
            list_small_polygons.append(polrow)

# split polygon with terminus line
list_areas = []
list_small_polygons = []
# Cryohackathon terminus lines
for (day, terminus_line) in zip(df_terminus.Date, df_terminus.geometry):
    newrow = {}
    newrow["Date"] = day.date()
    append_positions_by_section(newrow, terminus_line, list_small_polygons)
    list_areas.append(newrow)
# save polygons with glacier part of sections
gdf_small_polygons = gpd.GeoDataFrame(list_small_polygons, crs=df_sec.crs)
gdf_small_polygons.to_file(data_path+"section_polygons_glacier_only.gpkg", driver="GPKG")

# add McNabb terminus lines (for long term trend, don't help much for seasonal variability)
# for (yr, day, terminus_line) in zip(df_terminus_mcnabb.year, df_terminus_mcnabb.doy, df_terminus_mcnabb.geometry):
#     newrow = {}
#     day = str(int(day)+1)              # to handle days = 0 which exist in this dataset
#     newrow["Date"] = dt.datetime.strptime(yr + "-" + day, "%Y-%j").date()
#     append_positions_by_section(newrow, terminus_line)
#     list_areas.append(newrow)

# save all relative terminus positions
df_areas = pd.DataFrame(list_areas)
df_areas = df_areas.sort_values(by=["Date"])

# prepare lowess filter to smooth time series, easier on day index rather than datetime format
t0 = dt.datetime.strptime("2017-01-01", "%Y-%m-%d").date()
t_index_data   = np.zeros(len(df_areas.Date), dtype=int)
for (i, da) in enumerate(df_areas.Date):
    deltat = da - t0
    t_index_data[i] = deltat.days

t_smooth = pd.date_range(start=df_areas.Date[0], end=df_areas.Date[len(df_areas.Date)-1], freq="10D").date
t_index_smooth = np.zeros(len(t_smooth), dtype=int)
for (i, da) in enumerate(t_smooth):
    deltat = da - t0
    t_index_smooth[i] = deltat.days

# smooth and plot
lowess = sm.nonparametric.lowess
df_terminus_smooth = pd.DataFrame({"Date":t_smooth})
cols = ["brown", "pink", "darkgreen", "lightgreen", "grey"]    # only give four colours to only plot the first four sections
for (s,c) in zip(df_sec.section,cols):
    z = lowess(df_areas[s] - np.mean(df_areas[s]), t_index_data, frac=1/20, xvals=t_index_smooth)
    df_terminus_smooth[s+" [m]"] = z

# save smoothed time series
df_terminus_smooth.to_csv(data_path+"terminus_position_smooth.csv", index=False)
