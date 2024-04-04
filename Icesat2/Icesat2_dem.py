from sliderule import sliderule, icesat2
import geoutils as gu
import xdem
import argparse

# Step 2: Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Rasterize Icesat-2 Data")

# Step 3: Add an argument for 'res'
parser.add_argument('--res', type=float, help="Rasterization resolution")

# Step 4: Parse the command-line arguments
args = parser.parse_args()

# Step 5: Access the value of 'res' and use it in your script
if args.res is not None:
    res = args.res
    print(f"Value of 'res' provided from the command line: {res}")
else:
    # If 'res' is not provided as an argument, you can set a default value or handle it accordingly
    res = 50  # Default value

# Now, you can use the 'res' variable in your script as needed
print(f"Final value of 'res' in your script: {res}")


icesat2.init("slideruleearth.io", verbose=False)
mala = sliderule.toregion('/data/Sittlein/Coreg/AOI.geojson')

parms = icesat2.init("slideruleearth.io")
parms = {
    "srt": 0,
    "len": 40,
    "res": 20,
    "cnf": 4,
    "atl08_class": [],
    "maxi": 1,
    "ats": 20.0,
    "cnt": 10,
    "H_min_win": 3.0,
    "sigma_r_max": 5.0,
    "poly": [
        {
            "lon": -141.311646,
            "lat": 59.677726
        },
        {
            "lon": -139.680176,
            "lat": 59.677726
        },
        {
            "lon": -139.680176,
            "lat": 60.591283
        },
        {
            "lon": -141.311646,
            "lat": 60.591283
        },
        {
            "lon": -141.311646,
            "lat": 59.677726
        }
    ],
    "asset": "icesat2"
}
gdf = icesat2.atl06p(parms, asset="icesat2")


import pandas as pd
import numpy as np

# Convert the datetime index to days
gdf.set_index(pd.to_datetime(gdf.index).to_period('D'), inplace=True)

print(gdf)

# Group the DataFrame by the 'time' column
grouped = gdf.groupby('time')[['h_mean', 'geometry']]


from pathlib import Path
import os


# Define the output directory
savepath = Path('/data/Sittlein/Coreg/Icesat2')

# Create the output directory if it doesn't exist
os.makedirs(savepath, exist_ok=True)

# Reproject the geometry column to Alaska Albers
import pyproj
from shapely.ops import transform
wgs84 = pyproj.CRS('EPSG:4326')
utm = pyproj.CRS('EPSG:3338')
project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform

# Open reference DEM for the transform
ref_dem = xdem.DEM("/data/Sittlein/Coreg/Cop30.tif").reproject(dst_res = res)
ref_dem.set_nodata(np.nan)
# Get its transform
xform = ref_dem.transform

# Inverse the transform
xinv = ~xform



# Iterate over each group
for group, group_df in grouped:
    
    if len(group_df)>1:
        
        # Reproject the geometry column to Alaska Albers
        group_df.geometry = [transform(project, group_df.iloc[i].geometry) for i in range(len(group_df))]
        
        # Prepare the grids
        test = np.full(ref_dem.shape, 0)
        count = np.zeros(test.shape)
        
        # Populate the grid
        for i in range(len(group_df)):
            
            # Calculate the grid cell indices for the point
            idx = tuple(map(int, np.round(xinv*(group_df.iloc[i].geometry.x, group_df.iloc[i].geometry.y))))
            try:
                # Populate the grid by adding the points
                test[idx[1], idx[0]] += group_df.iloc[i].h_mean
                # Count the number of points in each cell
                count[idx[1], idx[0]] += 1
            except:
                continue
        # Divide the grid by the number of points in each cell
        count[count==0] = np.nan
        test = test/count
        

        # Define the output GeoTIFF file path
        output_tiff_path = f'{savepath}/{str(group_df.iloc[0].name)}_0.tif'


        gu.Raster.from_array(data=test, transform=ref_dem.transform, crs=ref_dem.crs, nodata=np.nan).save(output_tiff_path)
        

        print(f'DEM {str(group_df.iloc[0].name)} finished')
        
    else:
        print(f'DEM {str(group_df.iloc[0].name)} Failed (only 1 point)')
