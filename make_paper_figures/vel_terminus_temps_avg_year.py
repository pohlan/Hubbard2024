import os
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import statsmodels.api as sm
from pathlib import Path

# define data paths relative to this script's location
script_dir = Path(__file__).resolve().parent
project_root = script_dir.parent
terminus_data_path = project_root / "data" / "terminus_data" / ""  # folder where all the data for terminus position and terminus sections is stored
veloc_data_path    = project_root / "data" / "velocity" / ""
weather_data_path  = project_root / "data" / "weather" / ""
mooring_data_path  = project_root / "data" / "mooring" / "GAK" / ""
seismic_data_path  = project_root / "data" / "seismic" / ""
figure_path        = project_root / "make_paper_figures" / ""

# load datasets
df_veloc      = pd.read_csv(veloc_data_path / "velocity_by_section.csv")
df_terminus   = pd.read_csv(terminus_data_path / "retreat_by_section.csv")
df_ablation   = pd.read_csv(terminus_data_path / "frontal_ablation.csv")
df_retr_rate  = pd.read_csv(terminus_data_path / "retreat_rate_by_section.csv")
df_weather    = pd.read_csv(weather_data_path / "compiled_weather_data_PDDs.csv")
df_mooring    = pd.read_csv(mooring_data_path / "GAK_smooth.csv")
df_seismic    = pd.read_csv(seismic_data_path / "hubval_gq_1988-2024.csv", names=['lat', 'long', 'depth', 'Date', 'evid', 'mi'], header=1)

# bring time vectors into right format
t                  = [dt.datetime.strptime(day, "%Y-%m-%d").date() for day in df_terminus.Date]
df_veloc.Date      = pd.to_datetime(df_veloc.Date)
df_terminus.Date   = pd.to_datetime(df_terminus.Date)
df_ablation.Date   = pd.to_datetime(df_ablation.Date)
df_retr_rate.Date   = pd.to_datetime(df_retr_rate.Date)
df_mooring.Date    = pd.to_datetime(df_mooring.Date)
df_weather.rename(columns={"Unnamed: 0":"Date"}, inplace=True)
df_weather.Date = pd.to_datetime(df_weather.Date)  # already has a sample each day
df_seismic.Date = pd.to_datetime(df_seismic.Date)

# resample/interpolate time series to daily values (so they can later be averaged over the year more easily) and average over section
df_veloc = df_veloc.resample("d", on="Date").mean().interpolate()
df_terminus = df_terminus.resample("d", on="Date").mean().interpolate()
df_ablation = df_ablation.resample("d", on="Date").mean().interpolate()
df_retr_rate = df_retr_rate.resample("d", on="Date").mean().interpolate()
df_mooring = df_mooring.resample("d", on="Date").mean().interpolate()
df_seismic = df_seismic.resample("d", on="Date").count()
df_weather.set_index("Date",inplace=True)

# use lowess filter to smooth temperatue time series, easier on day index rather than datetime format
t0 = pd.Timestamp(dt.datetime.strptime("2017-01-01", "%Y-%m-%d").date())
t_index_data   = np.zeros(len(df_weather.index), dtype=int)
for (i, da) in enumerate(df_weather.index):
    deltat = da - t0
    t_index_data[i] = deltat.days
lowess = sm.nonparametric.lowess
z = lowess(df_weather["Yakutat_AirTemp_C"], t_index_data, frac=1/200, xvals=t_index_data)
df_weather["AirTempC_Yakutat_smooth"] = z

# # plot temperature
# plt.figure(figsize=(20,5))
# plt.plot(df_weather.index, df_weather["Yakutat_AirTemp_C"], label="Yakutat original")
# plt.plot(df_weather.index, df_weather["AirTempC_Yakutat_smooth"], lw=3, label="Yakutat smooth")
# plt.legend()

# define function to cut all the time series to the minimal overlap
def equal_period(df_inputs):
    dstart = np.max([np.min(df.index) for df in df_inputs])
    dend   = np.min([np.max(df.index) for df in df_inputs])
    dfs = [df.iloc[np.where((df.index >= dstart) & (df.index <= dend))[0]] for df in df_inputs]
    return dfs

# combine all data together in one dataframe
df_mooring1, df_veloc1, df_terminus1, df_ablation1, df_retr_rate1, df_weather1, df_seismic1    = equal_period([df_mooring, df_veloc, df_terminus, df_ablation, df_retr_rate, df_weather, df_seismic])

assert len(df_mooring1) == len(df_veloc1) == len(df_terminus1) == len(df_weather1) == len(df_seismic1)

df_all  = pd.DataFrame()
for (df, nam) in zip([df_veloc1, df_terminus1, df_ablation1, df_retr_rate1, df_weather1, df_mooring1.iloc[:,0:4], df_seismic1], ["velocity", "terminus", "ablation", "retreat_rate", "", "", "seis_events"]):
    if df_all.empty:
        df_all["Date"] = pd.to_datetime(df.index)  # df.index
    if (nam == "velocity") or (nam == "terminus") or (nam == "ablation") or (nam == "retreat_rate"):
        df_all[nam] = df.mean(axis=1).values
    else:
        df_all = df.merge(df_all, how="inner", on="Date")
df_save = df_all.rename(columns={"terminus":"terminus_retreat_m", "velocity":"surface_speed_m_per_yr", "ablation":"frontal_ablation_m_per_yr", "retreat_rate":"retreat_rate_m_per_yr", "Temperature [deg C] @ 25m":"OceanTemp@25m_degC", "Temperature [deg C] @ 50m":"OceanTemp@50m_degC", "Temperature [deg C] @ 100m":"OceanTemp@100m_degC", "Yakutat_Precip_mm":"Precip_Yakutat_24h_mm"})
df_save.to_csv("terminus_speed_weather_ocean.csv", index=False)

# resample precipitation data
dt_resample = 20
df_sum = df_save.resample(str(dt_resample)+"d", on="Date").mean().interpolate()

###########################
# 6-year timesieries plot #
###########################

fs = 20
plt.rc('font', size=fs)
fig, ax = plt.subplots(7,1, layout='constrained', figsize=(15,19))

# function for plotting velocity in blue in each panel
def plot_vel(axi):
    alpha = 0.8
    secax = axi.twinx()
    secax.plot(df_save["Date"], df_save.surface_speed_m_per_yr / 365, color="royalblue", lw=2, alpha=alpha)
    secax.set_ylabel("Surf. speed (m/d)", color="royalblue", alpha=alpha)
    secax.tick_params(axis='y', colors='royalblue')
    plt.setp(secax.get_yticklabels(), alpha=alpha)
    axi.xaxis.set_major_locator(mdates.YearLocator())
    axi.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

# terminus position
i = 0
ax[i].plot(df_save["Date"], df_save.terminus_retreat_m, color="black", lw=2)
ax[i].set_ylabel("Term. retreat (m)")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),410, 'a)')

# retreat rate
i = 1
ax[i].plot(df_save["Date"], df_save.retreat_rate_m_per_yr / 365, color="black", lw=2)
ax[i].set_ylim([-3,8])
ax[i].set_ylabel("Retreat rate (m/d)")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),6.5, 'b)')

# frontal ablation
i = 2
ax[i].plot(df_save["Date"], df_save.frontal_ablation_m_per_yr / 365, color="black", lw=2)
ax[i].set_ylabel("Front. abl. (m/d)")
ax[i].set_ylim([0,8])
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),6.5, 'c)')

# air temperature
i = 3
ax[i].plot(df_save["Date"], df_save.AirTempC_Yakutat_smooth, color="black", lw=2)
ax[i].set_ylabel("Air temp. (°C)")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),13, 'd)')

# precipitation
i = 4
ax[i].bar(df_sum.index, df_sum.Precip_Yakutat_24h_mm*1000, width=dt.timedelta(days=dt_resample), color="black", linewidth=0, alpha=0.3)
ax[i].set_ylabel("24hr precip. (mm)")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),900, 'e)')

# ocean temperature
i = 5
ax[i].plot(df_save["Date"], df_save["OceanTemp@50m_degC"], color="black")
ax[i].set_ylabel("Ocean temp. (°C)")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),10, 'f)')

# seismic events
i = 6
ax[i].plot(df_save["Date"], df_save["evid"], color="black")
ax[i].set_ylabel("Seismic events")
plot_vel(ax[i])
ax[i].text(dt.datetime(2016,12,15),9.5, 'g)')

# save
fig.savefig(figure_path / "six-years-timeseries.jpg")


#######################################################
# additional figures (Martin's Stanford presentation) #
#######################################################

# 1) plot velocity and terminus over one year (2017/2018)
fig, ax = plt.subplots(figsize=(15,4))
ax.plot(df_save["Date"], df_save.surface_speed_m_per_yr / 365, color="black", lw=2)
ax.set_ylabel("Speed (m/d)")
secax = ax.twinx()
secax.plot(df_save["Date"], df_save.terminus_retreat_m, color="royalblue", lw=2)
secax.set_ylabel("Terminus Retreat (m)", color='blue')
secax.set_yticks([200,300,400])
secax.set_yticklabels(labels=['200', '300', '400'],color='blue')
ax.set_xlim(pd.to_datetime('2017-08-01'), pd.to_datetime('2018-09-30'))
fig.savefig(figure_path / "vel_terminus_2017.jpg")

# 2) plot velocity and frontal ablation from 2017 to 2023, with green and red fill_between depending on which one is bigger
fig, ax = plt.subplots(figsize=(15,4))
ax.plot(df_save["Date"], df_save.surface_speed_m_per_yr / 365, color="black", lw=2)
ax.set_ylim([0.1,7.9])
ax.set_ylabel("Speed (m/d)")
secax = ax.twinx()
secax.plot(df_save["Date"], df_save.frontal_ablation_m_per_yr / 365, color="blue", lw=2)
secax.set_ylim([.1,7.9])
secax.set_ylabel("Frontal Abl. (m/d)", color='blue')
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.fill_between(df_save["Date"], df_save.surface_speed_m_per_yr / 365, df_save.frontal_ablation_m_per_yr / 365,
                df_save.surface_speed_m_per_yr>df_save.frontal_ablation_m_per_yr, color='green', alpha=0.2)
ax.fill_between(df_save["Date"], df_save.surface_speed_m_per_yr / 365, df_save.frontal_ablation_m_per_yr / 365,
                df_save.surface_speed_m_per_yr<df_save.frontal_ablation_m_per_yr, color='red', alpha=0.2)
fig.savefig(figure_path / "vel_frontal_ablation_2017_2023.jpg")


########################
# plot of average year #
########################

# make sure the end of the time series stops at the day of the year it starts (otherwise there are jumps in computing means over different years if days in e.g. January have less samples than in May)
end_year      = df_all.Date.iloc[-1].year
start_date    = df_all.Date.iloc[77].date()
end_date      = start_date.replace(year=end_year)
id_full_years = np.where((df_all.Date >= pd.Timestamp(start_date)) & (df_all.Date < pd.Timestamp(end_date)) & (df_all.Date != pd.Timestamp(2020, 2, 29)))  # removing leap day of 2020, makes sure there is the same amount of days in all of the years
df_full_years = df_all.iloc[id_full_years[0],:]
df_full_years.reset_index(inplace=True)

# check that there is the same amount of days for each day of the year
tind    = np.array([di.timetuple().tm_yday for di in df_full_years["Date"] ])
id_mv   = np.where((df_full_years.Date >= pd.Timestamp(2020,3,1)) & (df_full_years.Date <= pd.Timestamp(2020,12,31)))   # move all days of year after leap day one day forward
tind[id_mv] -= 1

for d in range(np.min(tind), np.max(tind)+1):
    assert np.sum(np.array(tind) == d) == end_year - start_date.year

def unique(x):
    return list(set(x))

# average over the year
ti_start_date = start_date.timetuple().tm_yday -1
strt_date = dt.date(start_date.year, 1, 1)
df_avyr = pd.DataFrame()
for col in df_full_years:
    if col == "Date":
        continue
    vls = []
    for dd in range(np.min(tind), np.max(tind)+1):
        val = np.mean(df_full_years[col][np.where(np.array(tind) == dd)[0]])
        vls.append(val)
    if df_avyr.empty:
        ds = []
        for doy in np.sort(unique(tind)):
            if doy <= ti_start_date:
                ds.append(strt_date.replace(year=strt_date.year+1)  + dt.timedelta(days=int(doy) - 1))
            else:
                ds.append(strt_date + dt.timedelta(days=int(doy) - 1))
        df_avyr["Date"] = ds
    df_avyr[col] = vls
df_avyr.sort_values(by="Date", inplace=True)
df_avyr.Date = pd.to_datetime(df_avyr.Date)   # otherwise Date is not in the correct format for .resample

# resample precipitation data
dt_resample = 10
df_sum = df_avyr.resample(str(dt_resample)+"d", on="Date").mean().interpolate()

# plot
fs = 20
plt.rc('font', size=fs)

def format_ax(ax, panel, plot_vel=True, draw_legend=False):
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[2,4,6,8,10,12]))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    if plot_vel:
        secax = ax.twinx()
        secax.plot(df_avyr["Date"], df_avyr["velocity"]/365, color="royalblue", linestyle="--", label="Mean 2017-2022")
        secax.set_ylabel(r"Surface speed ($\mathrm{m\,d^{-1}}$)", color="royalblue")
        secax.tick_params(axis='y', colors='royalblue')
        ax.set_zorder(secax.get_zorder()+1)
        ax.set_frame_on(False)                      # see https://matplotlib.org/stable/gallery/text_labels_and_annotations/date.html
        if draw_legend:
            secax.legend(loc="center left")
    ax.annotate(panel, (0.02,0.85), xycoords='axes fraction', fontweight='bold', fontsize=fs-3)

# can replace df_avyr by df_all to only average over sections and not year
fig, ax = plt.subplots(5,1, layout='constrained', figsize=(14,17))

# velocity
i = 0
alphas = [0.2, 0.35, 0.5, 0.7, 1.0]
for ii in range(1, 6):
    yr0 = pd.Timestamp(start_date) + pd.Timedelta(days=365*(ii-1))
    yrend = pd.Timestamp(start_date) + pd.Timedelta(days=365*ii)
    if yrend.year == 2020:                       # if there is a leap day, pd.Timedelta will take this into account, but we want to ignore it
        yrend = yrend + pd.Timedelta(days=1)
    id_yr = np.where((df_full_years.Date >= yr0) & (df_full_years.Date < yrend))[0]
    ax[i].plot(df_avyr.Date, df_full_years["velocity"][id_yr] / 365, color="royalblue", alpha=alphas[ii-1], label=str(yr0.year)+"-"+str(yr0.year+1))
    ax[i].set_ylabel(r"Surface speed ($\mathrm{m\,d^{-1}}$)")
    ax[i].legend()
format_ax(ax[i], "a", plot_vel=False)

# Terminus position
i = 1
ax[i].plot(df_avyr["Date"], df_avyr.terminus, color="black", lw=2); ax[i].set_ylabel("Terminus rel. retreat (m)")
format_ax(ax[i], "b", draw_legend=True)

# Air temperature
i = 2
ax[i].plot(df_avyr["Date"], df_avyr.AirTempC_Yakutat_smooth, color="black", lw=2); ax[i].set_ylabel("Air temperature (°C)")
format_ax(ax[i], "c")

# Precipitation
i = 3
ax[i].bar(df_sum.index, df_sum.Yakutat_Precip_mm *1000, width=dt.timedelta(days=dt_resample), color="grey", alpha=0.6, linewidth=0); ax[i].set_ylabel("24 hr precipitation (mm)")
format_ax(ax[i], "d")

# GAK temperatures
i = -1
stls = ["-", "--", "-.", "-.", "."]
for (n,col) in enumerate(df_all.filter(like="Temperature [deg C] @", axis=1)):
    ax[i].plot(df_avyr["Date"], df_avyr[col], stls[n], label=col.split(" ")[-1], color="black"); ax[i].set_ylabel("GAK temperature (°C)")
    ax[i].legend()
    ax[i].set_xlabel("Month of the year")
format_ax(ax[i], "e", False)

fig.savefig(figure_path / "four-panel-yearly.jpg")

