#!/bin/bash

#SBATCH --nodes=4
#SBATCH --time=16:00:00
#SBATCH --job-name="NWA12_wave_Ivan"
#SBATCH --output=NWA12_o.%j
#SBATCH --error=NWA12_e.%j
#SBATCH --qos=normal
#SBATCH --partition=batch
#SBATCH --clusters=c6
#SBATCH --account=bil-coastal-gfdl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qx7402@princeton.edu


# Avoid job errors because of filesystem synchronization delays
sync && sleep 1
export FI_VERBS_PREFER_XRC=0
export FI_CXI_RX_MATCH_MODE=hybrid
srun --ntasks=600 --cpus-per-task=1 --export=ALL ../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6

