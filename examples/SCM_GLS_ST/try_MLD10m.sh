#!/bin/bash

module load gcp
echo "Model started:  " `date`
current_folder=$(basename "$PWD")

NUM_TOT=20  ###overall loop number=y location
rm -rf 20*nc *.nc
###modify MOM_override file
y_loc=(300 200 150 130 110 90 70 50 30 10 0 -20 -40 -60 -80 -100 -120 -140 -200 -300)
file_name=('05_021.nc' '05_031.nc' '05_036.nc' '05_038.nc' '05_040.nc' '05_042.nc' '05_044.nc' '05_046.nc' '05_048.nc' '05_050.nc' '05_051.nc' '05_053.nc' '05_055.nc' '05_057.nc' '05_059.nc' '05_061.nc' '05_063.nc' '05_065.nc' '05_071.nc' '05_081.nc')
