import xarray as xr
import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
import csv
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Load velociy annual mosaic data
ds = xr.open_dataset("NSIDC-0776_RGI01A_2018_V02.0.nc")

vx = ds['vx'] # m/a
vy = ds['vy'] # m/a 

x = ds['x'] 
y = ds['y']

t = np.linspace(0, 20, 101) # years

# Create interpolators for vx and vy
vx_interp = RegularGridInterpolator((y, x), vx.values, bounds_error=False, fill_value=np.nan)
vy_interp = RegularGridInterpolator((y, x), vy.values, bounds_error=False, fill_value=np.nan)

# for now just find using a single point 
y0 = 248732.6
x0 = -3300764.2

vec0 = [x0, y0] 

# Check interpolated velocity at initial point
vx0 = vx_interp((y0, x0))
vy0 = vy_interp((y0, x0))
# print(f"vx at initial point: {vx0}")
# print(f"vy at initial point: {vy0}")

def vel(vec, t):
    x, y = vec
    vx_val = vx_interp((y, x))
    vy_val = vy_interp((y, x))
    return [vx_val, vy_val]

sol = odeint(vel, vec0, t)

# Remove rows with NaNs from sol (probably off glacier)
sol = sol[~np.isnan(sol).any(axis=1)]

# Save sol to CSV with columns x and y
df = pd.DataFrame({'x': sol[:, 0], 'y': sol[:, 1]})
df.to_csv('single_flowline.csv', index=False)

plt.plot(sol[:, 0], sol[:, 1], 'b-', label='Flowline')
plt.scatter([x0], [y0], color='red', label='Start')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.show()