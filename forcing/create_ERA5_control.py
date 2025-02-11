import xarray as xr
import xcdat as xc
import cftime
import datetime as dt
import numpy as np

FORCING_DIR = './'

# Load all ERA5 forcing data
clim_time = cftime.num2date(range(0,1460*6,6),'hours since 2021-01-01 00:00:00',calendar='365_day')
ds = xc.open_mfdataset(FORCING_DIR+'ERA5/raw/forcing_*.nc')

# Calculate climatology for every 6 hours
ds_clim = []
for date in clim_time:
    print(date)

    # For a certain day and time, get all years from 1971 to 2020
    tms = []
    for year in range(1971,2021):
        tms.append(cftime.DatetimeProlepticGregorian(year, date.month, date.day, date.hour, 0, 0, 0, has_year_zero=True))
    
    # Average over all values for that day and time
    ds_clim.append(ds.sel(time=tms).mean(dim='time').expand_dims(time=[date]))

# Concatenate all days
ds_clim = xr.concat(ds_clim,dim='time',data_vars='minimal')

# Open 2021 forcing data as base file (helps with keeping forcing file consistent with ACE)
ds = xr.open_dataset(FORCING_DIR+'ERA5/raw/forcing_2021.nc')

# Save all forcing inputs
ds['DSWRFtoa'] = (('time','latitude','longitude'), ds_clim['DSWRFtoa'].values)
ds['ocean_fraction'] = (('time','latitude','longitude'), ds_clim['ocean_fraction'].values)
ds['sea_ice_fraction'] = (('time','latitude','longitude'), ds_clim['sea_ice_fraction'].values)
ds['surface_temperature'] = (('time','latitude','longitude'), ds_clim['surface_temperature'].values)
ds['global_mean_co2'] = (('time'), np.full((ds_clim['time'].size),ds_clim['global_mean_co2'].values.mean()))

# Add December 31st to dataset
ds_ = ds.copy()
ds_ = ds_.isel(time=slice(-1,None))
ds_['time'] = np.array([dt.datetime(2020,12,31,18)])

ds_dec = xr.concat([ds_,ds],dim='time',data_vars='minimal')

# Save forcing
ds_dec.to_netcdf(FORCING_DIR+'ERA5/control/forcing.nc')
