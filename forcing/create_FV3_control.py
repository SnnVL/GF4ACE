import xarray as xr
import datetime as dt

FORCING_DIR = './'

# Load forcing data
ds = xr.open_dataset('./fv3/forcing_2021.zarr')

# Add December 31st to dataset
ds_ = ds.copy()
ds_ = ds_.isel(time=slice(-1,None))
ds_['time'] = ds_['time'] - dt.timedelta(days=365)
ds_dec = xr.concat([ds_,ds],dim='time',data_vars='minimal')

# Save as NetCDF
ds_dec.to_netcdf('./fv3/control/forcing.nc')