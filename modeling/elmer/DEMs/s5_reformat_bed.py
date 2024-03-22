import numpy as np
import matplotlib.pyplot as plt
import rasterio
import pandas as pd

# open bed DEM 
file_name = '/Users/amyjenson/Documents/Github/Hubbard2024/modeling/elmer/DEMs/hubbard_bed.tif'
with rasterio.open(file_name) as src:
    zb = src.read(1)
    height = zb.shape[0]
    width = zb.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xd, yd = rasterio.transform.xy(src.transform, rows, cols)
    lons_bed= np.array(xd)
    lats_bed = np.array(yd)

zb = zb.flatten()
zb = np.nan_to_num(zb, nan=-9999)
value = -10000
mat = np.where(zb > value, True, False)
bedDEM = np.where(mat, zb, -9999)

xb = lons_bed.flatten()
yb = lats_bed.flatten()
bedDEM = bedDEM.flatten()

#Save the bed in an ascii file
#np.savetxt('hubbard_bed.dat',np.column_stack((xb, yb, bedDEM)), fmt="%10.4f %10.4f %10.4f")
np.savetxt('./hubbard_bed.dat',np.column_stack((xb, yb, bedDEM)),fmt="%10.4f %10.4f %10.4f")


#info for sif file under Grid2Dinterpololator for bed
x0 = np.min(xb)       #bottom left x coordinate
y0 = np.min(yb)       #bottom left y coordinate
lx = (np.max(xb) - np.min(xb))    #the x length of the covered domain
ly = (np.max(yb) - np.min(yb))    #the y length of the covered domain
nx = len(lons_bed) #lx/dx +1              # number of levels in x and y directions (n,m in the tables above)
ny = len(lats_bed) #ly/dy + 1               #the number of levels in x and y directions (n,m in the tables above)
print('bed, x0=', f'{x0:.4f}', 'y0=', f'{y0:.4f}', 'lx =', f'{lx:.4f}', 'ly=', f'{ly:.4f}','nx=', nx,'ny=', ny)
