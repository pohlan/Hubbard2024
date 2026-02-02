import geopandas as gpd
import numpy as np
import datetime as dt
import pandas as pd
import xarray
from shapely import Point

terminus_data_path = "/home/annegret/Projects/Hubbard2024/data/terminus_data/"        # folder where all the data for terminus position and terminus sections is stored
veloc_data_path    = "/home/annegret/Projects/Hubbard2024/data/velocity/"

def get_along_vector(geom):
    # poly = list(geom.geoms)  # make Polygons out of MultiPolygon to be able to extract coordinates
    x, y = geom.exterior.coords.xy
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
    i_edge = edge_length.index(max(edge_length))
    along_vector = np.array([x[i_edge+1]-x[i_edge], y[i_edge+1]-y[i_edge]])
    return along_vector

# load sections, terminus positions and velocity data
geodf_sections = gpd.read_file(terminus_data_path+"section_polygons_glacier_only.gpkg")
xds_veloc      = xarray.open_dataset(veloc_data_path+"velocity_Radar_2026.nc")
df_terminus    = pd.read_csv("data/terminus_data/terminus_position_smooth.csv")
# move terminus sections to same coordinates as velocity (the other way round is not as straightforward since vx and vy also depend on that)
geodf_sections.to_crs(xds_veloc.projection, inplace=True)

# v_t    = [dt.datetime.strptime(t.isoformat(), "%Y-%m-%dT%H:%M:%S") for t in xds_veloc["time"].values]  # time series of velocities

# velocity on same time steps as terminus positions
tnp_term   = [np.datetime64(d) for d in df_terminus["Date"]]
xds_interp = xds_veloc.interp(time=tnp_term)
xds_interp.rio.write_crs("epsg:"+xds_interp.projection, inplace=True)

# initialize dataframes
t_term     = [dt.datetime.strptime(d, "%Y-%m-%d").date() for d in df_terminus["Date"]]                     # time series of terminus positions
df_vel     = pd.DataFrame({"Date":t_term})
df_vel_ortho = pd.DataFrame({"Date":t_term})
df_advance = pd.DataFrame({"Date":t_term})
df_dLdt    = pd.DataFrame({"Date":t_term[1:-1]})
df_abl     = pd.DataFrame({"Date":t_term[1:-1]})
for (i,(sec,geom)) in enumerate(zip(geodf_sections.section, geodf_sections.geometry.values)):
    # get median velocity in the section
    clipped  = xds_interp.rio.clip([geom], geodf_sections.crs)     # extract velocity data points that lie within the section
    v        = np.sqrt(clipped.vx**2 + clipped.vy**2)                   # get absolute velocity
    v_median = v.median(dim=["x","y"])
    df_vel[sec+" median |v| [m/yr]"] = v_median.values         # save median velocity values in dictionary to save later

    # calculate frontal ablation
    along_vector = get_along_vector(geom)
    vx_median = clipped.vx.median(dim=["x","y"])
    vy_median = clipped.vy.median(dim=["x","y"])
    vel_flux  = np.zeros(len(vx_median))
    for (j,(vx, vy)) in enumerate(zip(vx_median, vy_median)):
        median_v_vector = [vx, vy]
        orthogonal_vector = (np.dot(along_vector, median_v_vector)/sum(along_vector**2))*along_vector # orthogonal to section box, roughly orthogonal to terminus line
        vel_flux[j] = np.sqrt(sum(orthogonal_vector**2))
    df_vel_ortho[sec+" orthog. |v| [m/yr]"] = vel_flux

    # terminus advance position
    advance = df_terminus[sec+" [m]"]-np.min(df_terminus[sec+" [m]"])
    df_advance[sec+" advance position [m]"] = advance

    # time derivative of advanced position -> advance rate in m/yr
    dts = np.diff(t_term)  # time interval of data
    assert all(dts == dts[0])
    dL_dt = 0.5*(np.diff(advance[0:-1]) + np.diff(advance[1:])) *365/dts[0].days
    # dL_dt = 0.5*(np.diff(df_terminus[sec+" [m]"][0:-1]) + np.diff(df_terminus[sec+" [m]"][1:])) *365/5 # derivative of advance position
    df_dLdt[sec+" advance rate [m/yr]"] = dL_dt

    # abl = v_median[1:-1] + dL_dt
    abl = vel_flux[1:-1] - dL_dt
    df_abl[sec+" frontal ablation [m/yr]"] = abl

df_vel.to_csv(veloc_data_path+"velocity_by_section.csv", index=False)
df_advance.to_csv(terminus_data_path+"advance_by_section.csv", index=False)
df_dLdt.to_csv(terminus_data_path+"advance_rate_by_section.csv", index=False)
df_abl.to_csv(terminus_data_path+"frontal_ablation.csv", index=False)
df_vel_ortho.to_csv(veloc_data_path+"velocity_ortho_by_section.csv", index=False)
