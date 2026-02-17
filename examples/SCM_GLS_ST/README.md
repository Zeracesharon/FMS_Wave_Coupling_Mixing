single column model case, use GLS (from GOTM package) but no Langmuir turbulence included. ST refers to shear-turbulence. To run this, you need build the Grid first using BuildExchangeGrid.csh, simply by ./BuildExchangeGrid.csh should work. For the settings of BuildExchangeGrid.csh, as we use simple cartesian grid, supergrid is used here, the overall domain is 3000m, the nlon and nlat is given by 8 as the NIGLOBAL NJGLOBAL for the SCM is 4. After obtaining the wind and wave forcing (at the folder of SCM_Forcing), you can then run the model. If you want to run once, use ./run.sh Or just run:  ../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6. If you want to run once in HPC, run: sbatch job_script.sh. Normally, SCM is very quick, you can just run it in login in nodes or interaction mode. If you want run multiple locations at once, use try.sh for 5mps, initial MLD=32m; try_10mps for 10mps, initial MLD=32m; try_MLD10m.sh for 5mps, initial MLD=10m;try_10mps_MLD10m.sh, initial MLD=10m; for which indicates different wind forcing files at different locations will be linked. If you only ocean_only driver which requires changing the location but not the wind forcing file, then use try_location.sh for both 5mps and 10mps. The corresponding MOM_override settings for idealized SCM hurricane cases are: 
10mps:
#override IDL_HURR_SCM_LOCY = -3.0E+05 (location changes)
#override IDL_HURR_TRAN_SPEED = 10.0
#override IDL_HURR_X0 = 1.296E+06
5mps:
#override IDL_HURR_X0 = 6.48E+05
#override IDL_HURR_TRAN_SPEED = 5.0
The command will be: nohup bash try.sh > submit_jobs.log 2>&1 & This will then run all of the cases and give you results folder indicates different location.

