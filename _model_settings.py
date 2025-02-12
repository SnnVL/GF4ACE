"""
Model settings for different versions of ACE.

Each model has the following settings:
- dir: directory name
- model_name: model name for plotting purposes
- ckpt_name: checkpoint file name
- lat_name: latitude variable name
- lon_name: longitude variable name
- tas_name: surface temperature variable name
- rsut_name: ToA shortwave radiation variable name
- rlut_name: ToA longwave radiation variable name
- rsdt_name: ToA incoming solar (SW) radiation variable name
- save_vars: list of variables to save from ACE run
- control_run_length: number of years for control run
- patch_run_length: number of years for patch run

"""

__author__ = "Senne Van Loon"
__date__ = "10 February 2025"

DIR = './'

model_dct = {
    'era': {
        'model_name' : 'ACE2-ERA5',
        'ckpt_name': 'ace2_era5_ckpt.tar',
        'lat_name' : 'latitude',
        'lon_name' : 'longitude',
        'tas_name': 'surface_temperature',
        'rsut_name': 'USWRFtoa',
        'rlut_name': 'ULWRFtoa',
        'rsdt_name': 'DSWRFtoa',
        #
        'save_vars': "[ULWRFtoa,USWRFtoa]",
        #
        'control_run_length': 20,
        'patch_run_length': 12,
    },
    'eam': {
        'model_name' : 'ACE-EAM',
        'ckpt_name': 'ema_ckpt.tar',
        'lat_name' : 'lat',
        'lon_name' : 'lon',
        'tas_name': 'TS',
        'rsut_name': 'top_of_atmos_upward_shortwave_flux',
        'rlut_name': 'FLUT',
        'rsdt_name': 'SOLIN',
        #
        'save_vars': "['FLUT','top_of_atmos_upward_shortwave_flux']",
        #
        'control_run_length': 20,
        'patch_run_length': 12,
    },
    'fv3': {
        'model_name' : 'ACE-FV3',
        'ckpt_name': 'ace_ckpt.tar',
        'lat_name' : 'latitude',
        'lon_name' : 'longitude',
        'tas_name': 'surface_temperature',
        'rsut_name': 'USWRFtoa',
        'rlut_name': 'ULWRFtoa',
        'rsdt_name': 'DSWRFtoa',
        #
        'save_vars': "['ULWRFtoa', 'USWRFtoa']",
        #
        'control_run_length': 20,
        'patch_run_length': 12,
    }
}