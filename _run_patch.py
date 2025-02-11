"""
Run patch simulations

USAGE:
python run_forcing.py <model> <patch_location> <patch_amplitude>

<model> = 'era', 'eam', or 'fv3' (or create new entry in _model_settings)
"""

__author__ = "Senne Van Loon"
__date__ = "10 February 2025"

import subprocess as sp
import os
import xarray as xr
import sys

# Import model settings
from _model_settings import DIR, model_dct

# Get model settings, patch location, and patch amplitude
print(sys.argv)
model = str(sys.argv[1])
mod_dct = model_dct[model]
patch_loc = str(sys.argv[2])
A = str(sys.argv[3])

##### Define directories

# Configuration directory
CONF_DIR = DIR+'configs/'+model+'/'

# Forcing directory
FORCING_DIR = DIR+'forcing/'+model+'/'+patch_loc+'/'
if not os.path.isdir(FORCING_DIR):
    sp.run("mkdir "+FORCING_DIR, shell=True)

# Make patch forcing if it does not exist
if os.path.isfile(FORCING_DIR+'forcing.nc') or patch_loc == 'control':
    # Get initial time
    with xr.open_dataset(FORCING_DIR+'forcing.nc') as ds:
        t0 = ds['time'][0]
        t0 = t0.expand_dims({'sample':1})
        t0 = t0.reset_coords()
else:
    ds_patch = xr.open_dataset(\
        DIR+'forcing/'+model+'/GFMIP_patches_'+A+'K.nc',\
        decode_times=False\
    )[patch_loc]
    with xr.open_dataset(FORCING_DIR+'../control/forcing.nc') as ds:
        # Save patch forcing
        ds['surface_temperature'] += ds_patch
        ds.to_netcdf(FORCING_DIR+'forcing.nc')

        # Get initial time
        t0 = ds['time'][0]
        t0 = t0.expand_dims({'sample':1})
        t0 = t0.reset_coords()

# Output directory
OUTPUT_DIR = DIR+'output/'+model+'/'+patch_loc+'/'
if not os.path.isdir(OUTPUT_DIR):
    sp.run("mkdir "+OUTPUT_DIR, shell=True)

# Initial condition directory
INIT_DIR = DIR+'initialization/'+model+'/'+patch_loc+'/'
if not os.path.isdir(INIT_DIR):
    sp.run("mkdir "+INIT_DIR, shell=True)

# Copy initial conditions
if patch_loc == 'control':
    sp.run("cp "+INIT_DIR+'../20191231_1800.nc '+INIT_DIR+'20191231_1800.nc', shell=True)

    # First year is used for spinup
    start_year = 2020
    run_length = mod_dct['control_run_length']+1
else:
    sp.run("cp "+INIT_DIR+'../control/20401231_1800.nc '+INIT_DIR+'20201231_1800.nc', shell=True)

    start_year = 2021
    run_length = mod_dct['patch_run_length']

##### Run simulations
for year in range(start_year,start_year+run_length):
    yr_str = str(year)
    if os.path.isfile(INIT_DIR+yr_str+'1231_1800.nc'):
        print("Skip "+yr_str)
        continue
    print("******* "+yr_str)

    # Make output path
    if not os.path.isdir(OUTPUT_DIR+yr_str):
        sp.run("mkdir "+OUTPUT_DIR+yr_str, shell=True)

    # Create configuration file
    with open(CONF_DIR+'gf_template.yaml','r') as f:
        conf = f.read()
        conf = conf.replace("%DIR",model)
        conf = conf.replace("%CKPT",mod_dct['ckpt_name'])
        conf = conf.replace("%PATCH",patch_loc)
        conf = conf.replace("%YEAR",yr_str)
        conf = conf.replace("%INITYEAR",str(year-1))
        conf = conf.replace("%SAVE_VARS",mod_dct['save_vars'])
    with open(CONF_DIR+patch_loc+'.yaml','w') as f:
        f.write(conf)

    # Run inference
    sp.run("python -m fme.ace.inference "+CONF_DIR+patch_loc+".yaml",shell=True,cwd=DIR)
    
    # Create new initialization
    ds_init = xr.open_dataset(OUTPUT_DIR+yr_str+'/restart.nc')
    ds_init['time'] = (('sample'),t0['time'].data)
    ds_init.to_netcdf(INIT_DIR+yr_str+'1231_1800.nc')

# Remove forcing
sp.run("rm -r "+FORCING_DIR,shell=True)