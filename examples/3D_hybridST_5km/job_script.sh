#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --job-name="3D_HURR_GLS25km"
#SBATCH --output=3D_HURR_o.%j
#SBATCH --error=3D_HURR_e.%j
#SBATCH --qos=normal
#SBATCH --partition=batch
#SBATCH --clusters=c5
#SBATCH --account=gfdl_o
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qx7402@princeton.edu


# Avoid job errors because of filesystem synchronization delays
sync && sleep 1
export FI_VERBS_PREFER_XRC=0
srun --ntasks=100 --export=ALL ../../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6
