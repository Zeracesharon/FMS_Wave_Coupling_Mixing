This example provides the matched 1D experiments (corresponding to the 3D cases with the same translation speed). The workflow is basically the same as the 3D_hybrid_ST_5km case, but the grid build must use a single-column model (SCM), and thus the wind-forcing input files will be linked to winds at multiple locations.

You will need to run BuildExchangeGrid.csh first to build the ocean grid, and then check the settings in MOM_override and SIS_override (especially the NIGLOBAL and NJGLOBAL values), as well as the hybrid-scheme settings. Then use srun.sh or srun_20.sh to run the job.

The difference between srun.sh and srun_20.sh is the number of locations (srun.sh: 20 points; srun_20.sh: 360 points) and the corresponding wind-forcing file name.
