import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.ops import split
import datetime as dt

# read in terminus traces (UTM07)
terminus_path = 'terminus_data/Hubbard_2017_2021_UTM07.shp'
df_terminus = gpd.read_file(terminus_path)   # df_terminus.geometry is already in LineString format

# read in section polygons
path_sections = 'terminus_data/merged_sections.shp'
df_sec = gpd.read_file(path_sections) # df_sA.geometry is already in Polygon format
df_sec = df_sec.to_crs(df_terminus.crs)

# split polygon with terminus line
list_areas = []
n = 0
for (day, terminus_line) in zip(df_terminus.Date, df_terminus.geometry):
    newrow = {}
    newrow["Date"] = dt.datetime.strptime(day, "%Y-%m-%d").date()
    for (sect_name, sect_geom) in zip(df_sec.layer, df_sec.geometry):
        split_pols = split(sect_geom, terminus_line)
        if len(list(split_pols.geoms)) != 2:
            # print(n)
            continue
        assert list(split_pols.geoms)[0].area > list(split_pols.geoms)[1].area
        newrow[sect_name] = list(split_pols.geoms)[0].area
    list_areas.append(newrow)
    n += 1
df_areas = pd.DataFrame(list_areas)
df_areas = df_areas.sort_values(by=["Date"])

plt.plot(df_areas.Date,df_areas.section_A / np.max(df_areas.section_A), label="section A")
plt.plot(df_areas.Date,df_areas.sectionB_2 / np.max(df_areas.sectionB_2), label="section B2")
plt.plot(df_areas.Date,df_areas.sectionD / np.max(df_areas.sectionD), label="section D")
plt.legend()
plt.savefig("terminus_vs_time.jpg")
