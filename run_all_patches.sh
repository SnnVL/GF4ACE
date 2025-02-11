#!/bin/bash

#### File to run all patches for a given model and patch amplitude ####
# USAGE: ./run_all_patches.sh <model> <patch_amplitude>


# Model & patch amplitude
model=$1         # "era", "fv3", or "eam"
A=$2             # Amplitude of patch in K

# Create patches if not already created
if [ ! -f "./forcing/${model}/GFMIP_patches_${A}K.nc" ]; then
    python -u create_patches.py ${model} ${A}
fi

# All patches
declare -a arr=("patch_180.0E-80.0N" "patch_50.35E-80.0N" "patch_160.0E-70.0N" "patch_277.0E-70.0N" "patch_33.9E-70.0N" "patch_180.0E-60.0N" "patch_260.0E-60.0N" "patch_340.0E-60.0N" "patch_60.0E-60.0N" "patch_140.0E-60.0N" "patch_160.0E-50.0N" "patch_222.2E-50.0N" "patch_284.5E-50.0N" "patch_346.7E-50.0N" "patch_48.92E-50.0N" "patch_111.1E-50.0N" "patch_180.0E-40.0N" "patch_232.2E-40.0N" "patch_284.4E-40.0N" "patch_336.6E-40.0N" "patch_28.87E-40.0N" "patch_81.08E-40.0N" "patch_133.3E-40.0N" "patch_160.0E-30.0N" "patch_200.0E-30.0N" "patch_240.0E-30.0N" "patch_280.0E-30.0N" "patch_320.0E-30.0N" "patch_0.0E-30.0N" "patch_40.0E-30.0N" "patch_80.0E-30.0N" "patch_120.0E-30.0N" "patch_180.0E-20.0N" "patch_220.0E-20.0N" "patch_260.0E-20.0N" "patch_300.0E-20.0N" "patch_340.0E-20.0N" "patch_20.0E-20.0N" "patch_60.0E-20.0N" "patch_100.0E-20.0N" "patch_140.0E-20.0N" "patch_160.0E-10.0N" "patch_200.0E-10.0N" "patch_240.0E-10.0N" "patch_280.0E-10.0N" "patch_320.0E-10.0N" "patch_0.0E-10.0N" "patch_40.0E-10.0N" "patch_80.0E-10.0N" "patch_120.0E-10.0N" "patch_180.0E0.0N" "patch_220.0E0.0N" "patch_260.0E0.0N" "patch_300.0E0.0N" "patch_340.0E0.0N" "patch_20.0E0.0N" "patch_60.0E0.0N" "patch_100.0E0.0N" "patch_140.0E0.0N" "patch_160.0E10.0N" "patch_200.0E10.0N" "patch_240.0E10.0N" "patch_280.0E10.0N" "patch_320.0E10.0N" "patch_0.0E10.0N" "patch_40.0E10.0N" "patch_80.0E10.0N" "patch_120.0E10.0N" "patch_180.0E20.0N" "patch_220.0E20.0N" "patch_260.0E20.0N" "patch_300.0E20.0N" "patch_340.0E20.0N" "patch_20.0E20.0N" "patch_60.0E20.0N" "patch_100.0E20.0N" "patch_140.0E20.0N" "patch_160.0E30.0N" "patch_200.0E30.0N" "patch_240.0E30.0N" "patch_280.0E30.0N" "patch_320.0E30.0N" "patch_0.0E30.0N" "patch_40.0E30.0N" "patch_80.0E30.0N" "patch_120.0E30.0N" "patch_180.0E40.0N" "patch_232.2E40.0N" "patch_284.4E40.0N" "patch_336.6E40.0N" "patch_28.87E40.0N" "patch_81.08E40.0N" "patch_133.3E40.0N" "patch_160.0E50.0N" "patch_222.2E50.0N" "patch_284.5E50.0N" "patch_346.7E50.0N" "patch_48.92E50.0N" "patch_111.1E50.0N" "patch_180.0E60.0N" "patch_260.0E60.0N" "patch_340.0E60.0N" "patch_60.0E60.0N" "patch_140.0E60.0N" "patch_160.0E70.0N" "patch_277.0E70.0N" "patch_33.9E70.0N" "patch_180.0E80.0N" "patch_50.35E80.0N")

# First, run control simulation
python -u run_patch.py ${model} "control" None

# Run all patches
for i in "${arr[@]}"
do
    python -u run_patch.py ${model} "${i}_${A}K" ${A}
done

# Save net ToA radiation
python _get_R_values.py ${model} ${A}