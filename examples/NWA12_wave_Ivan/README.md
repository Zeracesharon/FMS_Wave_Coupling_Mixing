This Hurricane Ivan example uses a two-way coupled model. Only the files that differ from the NWA12_Ivan case are provided here. The restart files are the same as those in the NWA12_Ivan case.

This case setup is run on C6. A 600-core configuration is the most efficient setting and is recommended for wave components of 50×36 or 25×24. The workflow here is the same as in 3D_hybrid_twoway_5km. Since the ocean grid and topography files have been provided, you only need to go to INPUT/EXCHANGE_GRIDS/, run BuildExchangeGrids.csh in interactive mode, and then run python WaveMosaics.py.

Then go to WW3/PreProc/ and run:

python GenGrid.py
../../../../build/ncrc6.intel23/ww3_grid/REPRO/ww3_grid


Remember to change line 38 of the ww3_grid.inp file to the correct NI and NJ values.

Finally, go to NWA12_wave_Ivan/ and modify ww3_multi.inp according to your settings, then run the model with:
nohup bash try_wave.sh > submit_jobs.log 2>&1 &
