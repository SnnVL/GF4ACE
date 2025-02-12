"""
File to run all patches simulations for a given model and patch amplitude

model = 'era', 'eam', or 'fv3'
A = patch amplitude in K

"""

__author__ = "Senne Van Loon"
__date__ = "11 February 2025"

model = 'ace'
A = 2

import xarray as xr
from module_GF import DIR, create_patches, run_patch, get_R_values

if __name__ == '__main__':
    # Create patches
    create_patches(model,A)

    # Run control
    run_patch(model,'control',A)

    # Get patch list
    with xr.open_dataset(\
            DIR+'forcing/'+model+'/GFMIP_patches_'+str(A)+'K.nc',decode_times=False\
        ) as ds_patches:
        patch_lst = [item for item in list(ds_patches.keys()) if item.startswith('patch')]

    # Run patches
    for patch_loc in patch_lst:
        run_patch(model,patch_loc,A)

    # Save ToA radiation
    get_R_values(model,A,patch_lst)