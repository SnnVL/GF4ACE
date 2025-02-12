# GF4ACE
Green's function experiments with the [Ai2 Climate Emulator (ACE)](https://github.com/ai2cm/ace)


## Usage

### Installation
To install the necessary packages, run
```
conda create --name fme python=3.10
conda activate fme
pip install fme
conda install -c conda-forge xcdat
pip install cartopy
pip install jupyterlab
conda install conda-forge::cmcrameri
```
See [ACE documentation](https://ai2-climate-emulator.readthedocs.io/en/latest/) for more information about ACE. 

### Data
Checkpoint and forcing data for ACE can be found at the [Ai2 Hugging Face repository](https://huggingface.co/collections/allenai/ace-67327d822f0f0d8e0e5e6ca4)

Initial conditions and output data for the Green's function simulations can be found on [Dryad](https://doi.org/10.5061/dryad.d2547d8cf)

### Running GF4ACE

First, download the necessary data for the model checkpoints, forcing files, and initial conditions. To create the climatology for ACE2-ERA5, run `./forcing/create_ERA5_control.py`. 

Forcing files should be stored in ```./forcing/<model>/control/forcing.nc``` 

Initial conditions should be stored in ```./initialization/<model>/20191231_1800.nc```

To get the net ToA radiation and to plot the Green's function, the incoming solar radiation and ocean masks should be stored in 
```
./output/<model>/rsdt.nc
./output/<model>/ocean_mask.nc
```

All model settings are stored in `_model_settings.py`. Make sure the checkpoint name is correct.

To perform all patch simulation for a certain `<model>` and `<patch_amplitude>`, change the model and amplitude in `_GF_driver.py` and run
```python _GF_driver.py```
Patch perturbations will be stored in `./forcing/<model>/GFMIP_patches_<patch_amplitude>K.nc` and the global-mean output stored in `./output/<model>/R_gm_monthly_<patch_amplitude>K.nc`. 

Finally, `make_GF_yearly.ipynb` contains the code to make the Green's functions and the reconstructions. For ease of comparison, the processed output for ACE2-ERA5, ACE-FV3, and ACE-EAM are available for download on [Dryad](https://doi.org/10.5061/dryad.d2547d8cf). To recreate the historical timeseries for $R$, we use the SST anomalies provided by [GFMIP](https://gfmip.org). The SST file should be stored in `./forcing/historical.nc`.


## Credits
This work is a collaborative effort between [Senne Van Loon](https://scholar.google.com/citations?user=6h7ft20AAAAJ&hl=en), [Maria Rugenstein](https://www.atmos.colostate.edu/people/faculty/rugenstein/), and [Elizabeth A. Barnes](https://barnes.atmos.colostate.edu). 

### References

Senne Van Loon, Maria Rugenstein, & Elizabeth A. Barnes (2025), Green's functions of ERA5 using the Ai2 Climate Emulator, preprint to be added

### License

This project is licensed under an MIT license.

MIT © [Senne Van Loon](https://github.com/SnnVL)
