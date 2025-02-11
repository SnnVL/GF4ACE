"""
Get net top-of-atmosphere radiation (R) values from ACE model output.

USAGE:
python _get_R_values.py <model> <patch_amplitude>

"""

__author__ = "Senne Van Loon"
__date__ = "10 February 2025"

import xarray as xr
import xcdat as xc
import sys
import cftime
import numpy as np

# Import model settings
from _model_settings import DIR, model_dct

# Get model settings, patch location, and patch amplitude
print(sys.argv)
model = str(sys.argv[1])
mod_dct = model_dct[model]
A = str(sys.argv[2])

mod_dct['vars'] = [mod_dct['rsut_name'],mod_dct['rlut_name']]

def get_global_avg(ds,var):

    # Weights
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = "weights"

    # Weighted mean
    da_weighted = ds[var].weighted(weights)
    global_mean = da_weighted.mean(("lon", "lat"), skipna=True)

    return global_mean

def get_monthly_data_ace(dct,patch,start_year,end_year):
    da = []
    for year in range(start_year,end_year+1):
        yr_str = str(year)
        
        # Load dataset
        ds = xc.open_dataset(\
            DIR+'output/'+model+'/'+patch+'/'+yr_str+'/monthly_mean_predictions.nc',\
            decode_times=False,\
            add_bounds=False\
        )
        ds = ds.isel(sample=0,time=slice(-2))
        ds.lat.attrs["units"] = "degrees_north"

        # Encode time coordinates
        ds['time'] = cftime.num2date([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],'days since '+yr_str+'-01-01 00:00:00',calendar='365_day')
        ds = ds.reset_coords()
        ds = ds.drop_vars('valid_time')
        ds = ds.drop_vars('init_time')

        # Add time bounds
        ds['time'].encoding['calendar'] = '365_day'
        ds = ds.bounds.add_missing_bounds()

        # Get ToA values
        ds = ds[dct['vars']]
        da.append(ds)
    return xr.concat(da,dim='time')

def get_R_ace(dct,patch,start_year=2021,end_year=2032):

    # Get all ToA data
    ds = get_monthly_data_ace(dct,patch,start_year,end_year)
    
    # Get incoming solar
    dy = end_year-start_year+1
    rsdt = xc.open_dataset(DIR+'output/'+dct['name']+'/rsdt.nc')
    ds[dct['rsdt_name']] = (('time','lat','lon'),np.tile(rsdt[dct['rsdt_name']].values,(dy,1,1)))

    # Total radiation
    ds['R_'+patch] = ds[dct['rsdt_name']] - ds[dct['rsut_name']] - ds[dct['rlut_name']]
    ds = ds['R_'+patch].to_dataset()

    ds_gm = get_global_avg(ds,'R_'+patch).to_dataset()

    return ds, ds_gm


# Get control
ds, ds_gm = get_R_ace(\
    mod_dct, 'control',\
    start_year=2021, end_year=2021+mod_dct['control_run_length']\
)

# Rename control time
ds = ds.rename({'time':'time_control'})
ds_gm = ds_gm.rename({'time':'time_control'})


# Load patch list
with xc.open_dataset(\
        DIR+'forcing/'+model+'/GFMIP_patches_'+A+'K.nc',decode_times=False\
    ) as ds_patches:
    patch_lst = [item for item in list(ds_patches.keys()) if item.startswith('patch')]

# Loop over all patches
for ii, patch_loc in enumerate(patch_lst):
    print(patch_loc)

    patch, patch_gm = get_R_ace(\
        mod_dct, patch_loc, \
        start_year=2021, end_year=2021+mod_dct['patch_run_length']\
    )
    ds = xr.merge([ds,patch])
    ds_gm = xr.merge([ds_gm,patch_gm])

# Save to file
ds.to_netcdf(DIR+'output/'+model+'/R_monthly_'+A+'K.nc')
ds_gm.to_netcdf(DIR+'output/'+model+'/R_gm_monthly_'+A+'K.nc')