This example is used for running a two-way coupled ocean–wave model for 3D idealized hurricane cases. Only the files that differ are provided in this folder; all other files are the same as in 3D_hybridST_5km.

We first need to set up the grids for the ocean, wave, and ocean–wave exchange. The steps are as follows:

In the case folder, i.e., 3D_hybrid_twoway_5km/, run ./BuildExchangeGrid_new.csh in interactive mode (because you need srun). You will see the ocean grid files in the 3D_hybrid_twoway_5km/INPUT/ directory. Since we choose the wave hgrid and topog.nc to be identical to the ocean grid, we directly use the ocean grid files for the wave component.

Then go to 3D_hybrid_twoway_5km/INPUT/EXCHANGE_GRIDS/ to set up the exchange grids between the ocean and wave. Using the existing symlinks, run ./BuildExchangeGrids.csh. You will obtain the exchange grids between different components, including the wave and ocean. Then run python WaveMosaics.py to obtain grid_spec_waves.nc. At this point, you have built the ocean grid, wave grid, and exchange grid.

Next, you need to update the wave grid under the WW3/PreProc/ folder. Go to 3D_hybrid_twoway_5km/WW3/PreProc/ and run python GenGrid.py. Make sure the grid handle and grid name in GenGrid.py are set correctly for your files. Then set up your ww3_grid.inp file, specifically line 38 for NI and NJ. After modifying this file, run:

../../../../build/ncrc6.intel23/ww3_grid/REPRO/ww3_grid


Now return to the folder 3D_hybrid_twoway_5km/, and you are ready to run the model. Similar to 3D_hybridST_5km, if you want to submit the job and send the output folder to the archive analysis location, use srun.sh; if you want to run in interactive mode, use single.sh.

Note that, compared to 3D_hybridST_5km, this job script includes additional post-processing for wave output files.
