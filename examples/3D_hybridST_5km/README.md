This example simulates 3D idealized cases using the hybrid ST scheme (excluding LT effects, so no wave information is required). The first step is to build the ocean grid using BuildExchangeGrid_new.csh. The current resolution is 5 km. The domain spans from −13.5° to 13.5° in x, and from −8.1° to 8.1° in y, corresponding to roughly 3000 km in x and 1800 km in y.

Using the supergrid, nx = 2 × (3000 km / 5 km) = 1200 and ny = 2 × (1800 km / 5 km) = 720. In this case, dx = dy = (5 km) / 2 = 2500 m. The bottom depth is 4000 m.

Please pay attention to the srun task count (-n): not all core counts work well. 30 works in our setup, and you can try other values. If it fails, change the number of tasks and retry until it succeeds.

Run the command: srun -n 30 ./BuildExchangeGrid_new.csh. You will obtain the grid files inside the ./INPUT/ directory.

Next, link the wind file into ./INPUT/ from your 3D wind-forcing files. Remember to make the symlink name consistent with the entry in data_table.

Then modify NI and NJ in MOM_input or MOM_override, and also in SIS_input/SIS_override:

#override NIGLOBAL = 600
#override NJGLOBAL = 360


If you want to build a different resolution, change BuildExchangeGrid_new.csh and also update MOM_override and SIS_override accordingly. There are examples of BuildExchangeGrid_new.csh to create grids at 25 km, 10 km, 1 km, and 500 m resolution. Whenever you change resolution, remember to update the NIGLOBAL and NJGLOBAL values in MOM_override and SIS_override.

Now that you have built the grid and linked the input wind file, you can run the model.

If you want to run in interactive mode, use: nohup bash single.sh > submit_jobs.log 2>&1 &. Modify line 2 of single.sh to change the name of the output folder.

If you want to submit a batch job, use: nohup bash srun.sh > submit_jobs.log 2>&1 &, which will submit the job defined in job_script.sh. Modify lines 5 and 44 to change the output folder name and the archive destination (for the analysis node), respectively.

Note that for different hurricane translation speeds, you only need to link the corresponding wind-forcing files (and wave-forcing files, if applicable) into the INPUT/ directory and add them to data_table.
