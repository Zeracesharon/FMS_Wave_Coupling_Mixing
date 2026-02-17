This example uses the hybrid approach to simulate Hurricane Ivan (2004) without LT effects and with a one-way coupling setup. The 1/12° Northwest Atlantic regional grid is used. The grid files, initial conditions, and boundary conditions are all prepared following Ross et al. (2023). The restart files are described in the hybrid-approach paper. Most preparations are already complete; you only need to run:
nohup sbatch try_uncou.sh > submit_jobs.log 2>&1 &.

Please modify job_script.sh to select the desired nodes and submission settings. If you change the core number, you will need to use different MOM_layout and SIS_layout; in that case, update try_uncou.sh lines 8 and 9. Modify lines 243 and 554 to set the destination paths for the output and results on the analysis node.

If you want to run one-way coupling for the Ivan model, simply replace INPUT/MOM_override with MOM_override_oneway (first rename it to MOM_override, then place it in the INPUT/ directory).

To correctly link all the grid files and initial conditions for the Ivan cases, you will need to create the following symlinks in the directory FMS_Wave_Coupling_Mixing/examples:

ln -sf /gpfs/f5/cefi/world-shared/datasets datasets
ln -sf /gpfs/f5/gfdl_med/world-shared/northwest_atlantic/nwa12_input datasets1
ln -sf /gpfs/f5/gfdl_med/world-shared/northwest_atlantic/era5 era5
ln -sf /gpfs/f5/gfdl_o/world-shared/oneway_Stokes/Stokes stokes


Reference
Ross, A.C., Stock, C.A., Adcroft, A., Curchitser, E., Hallberg, R., Harrison, M.J., Hedstrom, K., Zadeh, N., Alexander, M., Chen, W. and Drenkard, E.J., 2023. A high-resolution physical-biogeochemical model for marine resource applications in the northwest Atlantic (MOM6-COBALT-NWA12 v1.0). Geoscientific Model Development Discussions, 2023, pp.1–65.
