"""
Create patch perturbations

USAGE:
python create_patches.py <model> <patch_amplitude>

<model> = 'era', 'eam', or 'fv3' (or create new entry in _model_settings)
"""

__author__ = "Senne Van Loon"
__date__ = "10 February 2025"

import xarray as xr
import numpy as np
import sys

# Import model settings
from _model_settings import DIR, model_dct

# Get model settings and patch amplitude
print(sys.argv)
model = str(sys.argv[1])
mod_dct = model_dct[model]
A = int(sys.argv[2])

# Patch perturbations
def patch_sst(A,patch_loc,patch_w,lat,lon):
    # Patch_loc: (lat_p,lon_p)
    # Patch_w: (lat_width,lon_width)

    # Define grid; use periodic boundary conditions
    nlon = lon.size
    lon = np.concatenate([lon-360,lon,lon+360])
    X, Y = np.meshgrid(lon,lat)
    
    # Zero out values outside patch
    zero_lat = np.logical_or(lat < patch_loc[0]-patch_w[0]/2 , patch_loc[0]+patch_w[0]/2 < lat)
    zero_lon = np.logical_or(lon < patch_loc[1]-patch_w[1]/2 , patch_loc[1]+patch_w[1]/2 < lon)

    # Create patch perturbation
    sst = A*np.cos(np.pi/2*(X-patch_loc[1])/(patch_w[1]/2))**2 \
            *np.cos(np.pi/2*(Y-patch_loc[0])/(patch_w[0]/2))**2

    # Zero out values outside patch
    sst[zero_lat,:] = 0.
    sst[:,zero_lon] = 0.

    # Sum over periodic boundary conditions
    sst = sst[:,:nlon]+sst[:,nlon:2*nlon]+sst[:,-nlon:]

    return sst

# Patch label
def get_patch_str(lon_loc,lat_loc,A):

    def round_to_n(x,n):
        if x == 0:
            return x
        else:
            return np.round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    def num_lab(x,n):
        return str(round_to_n(x,n))
    
    return 'patch_'+num_lab(lon_loc,4)+"E"+num_lab(lat_loc,4)+"N_"+str(A)+"K"

# Load latitude and longitude from control
with xr.open_dataset(DIR+'forcing/'+model+'/control/forcing.nc') as f:
    lat = f[mod_dct['lat_name']]
    lon = f[mod_dct['lat_name']]
    tpe = f[mod_dct['tas_name']].dtype

# Create patch dataset
patches = xr.Dataset(coords={
    mod_dct['lat_name']:lat.values,
    mod_dct['lat_name']:lon.values
})
patches[mod_dct['lat_name']] = lat
patches[mod_dct['lat_name']] = lon

# Loop over latitude
for lat_loc in np.arange(-80,81,10,dtype='float32'):

    # Patch width and longitude step
    patch_w = (20,80)
    lon_step = 40

    # Adjust patch width and longitude step for high latitudes
    if np.abs(lat_loc)>30:
        patch_w = (20,80/np.cos(np.deg2rad(lat_loc)))
        lon_step = 40/np.cos(np.deg2rad(lat_loc))

    # Alternate starting longitude
    if np.abs(lat_loc) in [0,20,40,60,80]:
        lon_start = 180
    else:
        lon_start = 160

    # Loop over longitude
    for lon_loc in np.arange(lon_start,lon_start+321,lon_step,dtype='float32'):
        lon_loc = lon_loc % 360
        patch_loc = (lat_loc,lon_loc)

        # Patch name
        patch_str = get_patch_str(lon_loc,lat_loc,A)
        print(patch_str)

        # Create patch
        sst = patch_sst(A,patch_loc,patch_w,lat.values,lon.values).astype(tpe)

        # Save patch
        patches[patch_str] = ((mod_dct['lat_name'],mod_dct['lat_name']),sst)

# Save all patches
patches.to_netcdf(DIR+'forcing/'+model+'/GFMIP_patches_'+str(A)+'K.nc')