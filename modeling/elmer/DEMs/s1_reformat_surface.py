import numpy as np
import matplotlib.pyplot as plt
import rasterio
import pandas as pd


# open surface DEM 
file_name = '/Users/amyjenson/Documents/Github/Hubbard2024/modeling/elmer/DEMs/hubbard_surface.tif'
with rasterio.open(file_name) as src:
    zs = src.read(1)
    height = zs.shape[0]
    width = zs.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xd, yd = rasterio.transform.xy(src.transform, rows, cols)
    lons_surf= np.array(xd)
    lats_surf = np.array(yd)

xs = np.flip(lons_surf.flatten())
ys = np.flip(lats_surf.flatten())
zs = np.flip(zs.flatten())
value = 0
mat = np.where(zs > value, True, False)
surfDEM = np.where(mat, zs, -9999)

# # plot surfDEM
# plt.figure()
# plt.xlim(0,100)
# plt.plot(z)
# plt.title('surface DEM')
# plt.show()

#Save the surface in an ascii file 
np.savetxt('./hubbard_surf.dat',np.column_stack((xs, ys, surfDEM)),fmt="%10.4f %10.4f %10.4f")
  
#info for sif file under Grid2Dinterpololator for surface
x0 = xs[0]       #bottom left x coordinate
y0 = ys[0]       #bottom left y coordinate
lx = (np.max(xs) - np.min(xs))    #the x length of the covered domain
ly = (np.max(ys) - np.min(ys))    #the y length of the covered domain
nx = len(lons_surf) #lx/dx +1              # number of levels in x and y directions (n,m in the tables above)
ny = len(lats_surf) #ly/dy + 1               #the number of levels in x and y directions (n,m in the tables above)
print('surface, x0=', f'{x0:.4f}', 'y0=', f'{y0:.4f}', 'lx =', f'{lx:.4f}', 'ly=', f'{ly:.4f}', 'nx=', nx,'ny=', ny)

#print(min(xb))
#print(min(yb))
