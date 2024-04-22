from rasterio.plot import show
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.enums import Resampling
from pyproj import Transformer
from rasterio.windows import Window
import pandas as pd
from pyproj import Proj, transform
import geopandas as gpd
from matplotlib_scalebar.scalebar import ScaleBar

fs = 24 #font size
        
# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))

########################
# add a flowline
points = pd.read_csv('centerline_points_100m.csv')
transformer = Transformer.from_crs("epsg:3413", "epsg:32607") # UTM 6N
x, y = transformer.transform(points.X.to_numpy(),points.Y.to_numpy())
ax.plot(x,y,'r')

########################
# coordinates of points for velocity plots 
points = pd.read_csv('centerline_points_3000m.csv')

#ax.scatter(points.X[1], points.Y[1] , alpha=1, c='blue', s=100)
#ax.scatter(points.X[3], points.Y[3] , alpha=1, c='orange', s=100)
#ax.scatter(points.X[5], points.Y[5] , alpha=1, c='black', s=100)

# plug in these coordinates to victors velocity datacube




########################

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig('map_inset.png')