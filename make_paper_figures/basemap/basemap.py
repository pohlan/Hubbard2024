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
import matplotlib.cm as cm
import matplotlib.colors as mcolors

fs = 24  # font size


def pct_clip(array, pct=[2, 98]):
    array_min, array_max = np.nanpercentile(array, pct[0]), np.nanpercentile(
        array, pct[1]
    )
    clip = (array - array_min) / (array_max - array_min)
    clip[clip > 1] = 1
    clip[clip < 0] = 0
    return clip


with rio.open("Planet_Hubbard.tif") as src:
    with rio.open(
        "RGB_Temp.tif",
        "w+",
        driver="GTiff",
        dtype=rio.float32,
        count=3,
        crs=src.crs,
        width=src.width,
        height=src.height,
        transform=src.transform,
    ) as dst:
        V = pct_clip(src.read(1))
        dst.write(V, 1)
        V = pct_clip(src.read(2))
        dst.write(V, 2)
        V = pct_clip(src.read(3))
        dst.write(V, 3)

# Create a larger figure
fig, ax = plt.subplots(figsize=(16, 16))

window = rio.windows.Window(3000, 3000, 6500, 6500)
with rio.open("RGB_Temp.tif") as src2:
    data = src2.read(window=window)
    transform = rio.windows.transform(window, src2.transform)
    show(data, transform=transform, ax=ax)


########################

sites = np.arange(1, 21, 1)

# Use a colormap to assign unique colors
cmap = cm.get_cmap("tab20", len(sites))

for idx, site in enumerate(sites):
    points = pd.read_csv(f"flowline_{site}.csv")
    transformer = Transformer.from_crs("epsg:3413", "epsg:32607")  # UTM 6N
    points_X, points_Y = transformer.transform(points.x.to_numpy(), points.y.to_numpy())
    color = cmap(idx)
    ax.plot(points_X, points_Y, label=f"Flowline {site}", color=color)

ax.legend()

# Get all text objects in the figure
text_objs = plt.gcf().findobj(plt.Text)

# Change font size for all text objects
font_size = fs  # Change this to the font size you desire
for text_obj in text_objs:
    text_obj.set_fontsize(font_size)

plt.savefig("map_with_flowlines.png")
