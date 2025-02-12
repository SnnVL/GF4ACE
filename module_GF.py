"""
Module containing functions to run Green's function experiments with ACE

Functions:
- create_patches(model, A): creates all patches for a given model and patch amplitude
- run_patch(model,patch_loc,A): runs a patch simulation for a given model and patch
- get_R_values(model,A,patch_lst): gets monthly net top-of-atmosphere radiation for a given model and patch amplitude
- setup_figure(projection,nCols=1,nRows=1,size=(15,15),mask=True): sets up a figure with a certain projection

"""

__author__ = "Senne Van Loon"
__date__ = "11 February 2025"

import xarray as xr
import xcdat as xc
import numpy as np
import subprocess as sp
import os
import cftime

# Import model settings
from _model_settings import DIR, model_dct

##########################################################################
#########################    CREATING PATCHES    #########################
##########################################################################

def patch_sst(A,patch_loc,patch_w,lat,lon):
    """
    Defines a patch perturbation in sea surface temperature

    # Patch_loc: (lat_p,lon_p)
    # Patch_w: (lat_width,lon_width)
    """

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

def get_patch_str(lon_loc,lat_loc,A):
    """
    Creates a label for a certain patch location and amplitude

    # Patch_loc: (lat_p,lon_p)
    # Patch_w: (lat_width,lon_width)
    """

    def round_to_n(x,n):
        if x == 0:
            return x
        else:
            return np.round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    def num_lab(x,n):
        return str(round_to_n(x,n))
    
    return 'patch_'+num_lab(lon_loc,4)+"E"+num_lab(lat_loc,4)+"N_"+str(A)+"K"

def create_patches(model, A):
    """
    Creates & saves all patches according to the GFMIP protocol, 
    for a given model and patch amplitude
    """

    # Get model settings
    mod_dct = model_dct[model]

    # Load latitude and longitude from control
    with xr.open_dataset(DIR+'forcing/'+model+'/control/forcing.nc') as f:
        lat = f[mod_dct['lat_name']]
        lon = f[mod_dct['lon_name']]
        tpe = f[mod_dct['tas_name']].dtype

    # Create patch dataset
    patches = xr.Dataset(coords={
        mod_dct['lat_name']:lat.values,
        mod_dct['lon_name']:lon.values
    })
    patches[mod_dct['lat_name']] = lat
    patches[mod_dct['lon_name']] = lon

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
            patches[patch_str] = ((mod_dct['lat_name'],mod_dct['lon_name']),sst)

    # Save all patches
    patches.to_netcdf(DIR+'forcing/'+model+'/GFMIP_patches_'+str(A)+'K.nc')

##########################################################################
########################    PATCH SIMULATIONS    #########################
##########################################################################

def run_patch(model,patch_loc,A):
    """
    Runs a patch simulation for a given model and patch
    """

    # Get model settings
    mod_dct = model_dct[model]

    # Configuration directory
    CONF_DIR = DIR+'configs/'+model+'/'
    if not os.path.isdir(CONF_DIR):
        sp.run("mkdir "+CONF_DIR, shell=True)

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
            DIR+'forcing/'+model+'/GFMIP_patches_'+str(A)+'K.nc',\
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
        with open(DIR+'configs/gf_template.yaml','r') as f:
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
    if not patch_loc == 'control':
        sp.run("rm -r "+FORCING_DIR,shell=True)

##########################################################################
########################     SAVING FUNCTIONS     ########################
##########################################################################

def write_line(f, x, y):
    f.write(str(x))
    f.write(" ")
    f.write(str(y))
    f.write("\n")

def get_global_avg(ds,var):
    """
    Get global average of a variable in dataset ds
    """

    # Weights
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = "weights"

    # Weighted mean
    da_weighted = ds[var].weighted(weights)
    global_mean = da_weighted.mean(("lon", "lat"), skipna=True)

    return global_mean


def get_ocean_avg(ds,var,mod_dct):
    """
    Get average over ice-free ocean for a variable in dataset ds
    """

    # Get ocean area
    weights = xr.DataArray(mod_dct['ocean_area'],
        coords={
            mod_dct['lat_name']: mod_dct['lat'],
            mod_dct['lon_name']: mod_dct['lon'],
        },name='weights')

    da_weighted = ds[var].weighted(weights)
    ocean_mean = da_weighted.mean((mod_dct['lon_name'], mod_dct['lat_name']), skipna=True).values

    return ocean_mean

def get_monthly_data_ace(dct,patch,start_year,end_year):
    """
    Get monthly data from the output of ACE
    """

    da = []
    for year in range(start_year,end_year+1):
        yr_str = str(year)
        
        # Load dataset
        ds = xc.open_dataset(\
            DIR+'output/'+dct['name']+'/'+patch+'/'+yr_str+'/monthly_mean_predictions.nc',\
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
    """
    Get monthly net top-of-atmosphere radiation from a single patch simulations
    """

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

def get_R_values(model,A,patch_lst):
    """
    Get monthly net top-of-atmosphere radiation from a single patch simulations for a certain model and patch amplitude
    """

    # Get model settings
    mod_dct = model_dct[model]
    mod_dct['name'] = model
    mod_dct['vars'] = [mod_dct['rsut_name'],mod_dct['rlut_name']]

    # Get control
    ds, ds_gm = get_R_ace(\
        mod_dct, 'control',\
        start_year=2021, end_year=2021+mod_dct['control_run_length']-1\
    )

    # Rename control time
    ds = ds.rename({'time':'time_control'})
    ds_gm = ds_gm.rename({'time':'time_control'})

    # Loop over all patches
    for patch_loc in patch_lst:
        print(patch_loc)

        patch, patch_gm = get_R_ace(\
            mod_dct, patch_loc, \
            start_year=2021, end_year=2021+mod_dct['patch_run_length']-1\
        )
        ds = xr.merge([ds,patch])
        ds_gm = xr.merge([ds_gm,patch_gm])

    # Save to file
    ds.to_netcdf(DIR+'output/'+model+'/R_monthly_'+str(A)+'K.nc')
    ds_gm.to_netcdf(DIR+'output/'+model+'/R_gm_monthly_'+str(A)+'K.nc')


##########################################################################
########################    PLOTTING FUNCTIONS    ########################
##########################################################################
    
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath

# Define map projections
data_crs = ccrs.PlateCarree()
proj_Global = ccrs.EqualEarth(central_longitude=200)
proj_North = ccrs.NorthPolarStereo()
proj_South = ccrs.NorthPolarStereo()
proj_Global_PC = ccrs.PlateCarree(central_longitude=210)

def setup_figure(projection,nCols=1,nRows=1,size=(15,15),mask=True):
    """
    Setup a figure with a certain projection
    """

    # Define map projection
    circ = False
    if projection.lower() == 'north':
        map_proj = proj_North
        circ = True
    elif projection.lower() == 'south':
        map_proj = proj_South
        circ = True
    elif projection.lower() == 'global':
        map_proj = proj_Global
    elif projection.lower() == 'global_pc':
        map_proj = proj_Global_PC
    else:
        raise NotImplementedError("Projection not implemented.")

    # Define circle
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # Plot land lines
    land_feature = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='50m',
        facecolor=(0.8,0.8,0.8),
        edgecolor='k',
        linewidth=.25,
    )
    
    # Create figure
    fig = plt.figure(figsize=size)
    if nCols==1 and nRows==1:
        ax = fig.add_subplot(1, 1, 1, projection=map_proj)

        ax.coastlines('50m', linewidth=0.8)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.gridlines(
            draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
        )

        if circ:
            ax.set_boundary(circle, transform=ax.transAxes)
        if mask:
            ax.add_feature(land_feature)
    elif nCols > 1 and nRows > 1:
        ax = np.empty((nCols,nRows),dtype=object)
        iF = 1
        for jj in range(nRows):
            for ii in range(nCols):
                ax[ii,jj] = fig.add_subplot(nRows, nCols, iF, projection=map_proj)

                ax[ii,jj].coastlines('50m', linewidth=0.8)
                ax[ii,jj].tick_params(axis='both', which='major', labelsize=10)
                ax[ii,jj].gridlines(
                    draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
                )
                if circ:
                    ax[ii,jj].set_boundary(circle, transform=ax.transAxes)
                if mask:
                    ax[ii,jj].add_feature(land_feature)
                iF += 1
    elif nCols > 1:
        ax = np.empty((nCols),dtype=object)
        for ii in range(nCols):
            ax[ii] = fig.add_subplot(1, nCols, ii+1, projection=map_proj)

            ax[ii].coastlines('50m', linewidth=0.8)
            ax[ii].tick_params(axis='both', which='major', labelsize=10)
            ax[ii].gridlines(
                draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
            )
            if circ:
                ax[ii].set_boundary(circle, transform=ax.transAxes)
            if mask:
                ax[ii].add_feature(land_feature)
    elif nRows > 1:
        ax = np.empty((nRows),dtype=object)
        for ii in range(nRows):
            ax[ii] = fig.add_subplot(nRows, 1, ii+1, projection=map_proj)

            ax[ii].coastlines('50m', linewidth=0.8)
            ax[ii].tick_params(axis='both', which='major', labelsize=10)
            ax[ii].gridlines(
                draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
            )
            if circ:
                ax[ii].set_boundary(circle, transform=ax.transAxes)
            if mask:
                ax[ii].add_feature(land_feature)
     
    return fig, ax