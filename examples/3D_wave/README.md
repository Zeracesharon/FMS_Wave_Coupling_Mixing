To run the wave model, you first need to set up the wave grid in the Grid/ folder. You will need wave_hgrid.nc and wave_topog.nc, which can be copied directly from the 3D MOM6 ocean cases after running BuildExchangeGrid. This means that you need to build the 3D ocean grid first (see 3D_hurricane/). Copy ocean_hgrid.nc and ocean_topog.nc, rename them as the wave components, and put them in the 3D_wave/Grid/ directory. Then run: python GenGrid.py.

Next, set up ww3_grid.inp. In particular, pay attention to the first line:
1.1 0.04118 25 24 0.5
to set your desired wave components and initial frequency. Also note the line containing 'CURV' T 'NONE', followed by 600 360. These two numbers indicate nx/2 and ny/2 in your ocean_hgrid.nc.

Also, at this line:

&OUTS USSP = 1, IUSSP = 9,
STK_WN = 0.01,0.03,0.05,0.08,0.13,0.2,0.35,0.5,0.7


you can choose the number of wavenumbers you want; here we chose 9. You can choose more, and/or use a wider range of values. Then run:
../../../../build/ncrc6.intel23/ww3_grid/REPRO/ww3_grid (the exact path depends on your build directory).

Wind forcing preparation (for WW3)

Next, you will need the wind file for the wave field. Go to the 3D_wave/Input/ directory. You can use the wind file that provides forcing for the 3D hurricane ocean simulations, but make sure the coordinate variable names are compatible: WW3 does not recognize Lon and Lat, but only lon, lat, or longitude, latitude. The code for generating wind input files for the wave field is in wavefield_windInput.ipynb.

After you obtain the wind forcing for the wave field (e.g., Wind_5km_10mi_wave.nc), set up ww3_prnc.inp. For the dimension names, variable names, and file name:
lon lat ; u10 v10 ; Wind_5km_10mi_wave.nc
then run: ../../../../build/ncrc6.intel23/ww3_prnc/REPRO/ww3_prnc. You will obtain a file named wind.ww3.

If you have an ice file, follow the same workflow, but change ww3_prnc.inp to use 'ICE' 'LL' T T.

Up to this point, the wave-grid and forcing inputs are ready.

Running WW3 (ww3_multi)

Now go to the folder 3D_wave/Wave/ and set up ww3_multi.inp. In particular:

Line 7: input file type

Line 8: initial and end time

Line 10: 19000101 000000 600 20230607 000000 — the third number is the output frequency (every 600 s here)

Make sure line 12 includes USP

Then you can run the wave model by submitting run.sh: sbatch run.sh. Once it finishes, you will find a file out_grd.ww3.

Post-processing (ww3_ounf)

Next, set up ww3_ounf.inp:

Line 7: initial and final dates for output, and the output frequency

Make sure line 22 includes USP

Line 40: choose the post-processing output file name

Then run: ../../../../build/ncrc6.intel23/ww3_ounf/REPRO/ww3_ounf. You will obtain ww3.201104_usp.nc, which will later be used by 3D1D_waveInput.ipynb.
