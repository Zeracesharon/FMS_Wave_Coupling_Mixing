#!/bin/bash

module load gcp
echo "Model started:  " `date`
current_folder=$(basename "$PWD")

NUM_TOT=20  ###overall loop number=y location
rm -rf 20*nc *.nc
###modify MOM_override file
y_loc=(300 200 150 130 110 90 70 50 30 10 0 -20 -40 -60 -80 -100 -120 -140 -200 -300)
file_name=('10_021.nc' '10_031.nc' '10_036.nc' '10_038.nc' '10_040.nc' '10_042.nc' '10_044.nc' '10_046.nc' '10_048.nc' '10_050.nc' '10_051.nc' '10_053.nc' '10_055.nc' '10_057.nc' '10_059.nc' '10_061.nc' '10_063.nc' '10_065.nc' '10_071.nc' '10_081.nc')
NUM_TOT=${#y_loc[@]}

