Wind–Wave Forcing Input Files: How to Generate (SCM + 3D Hurricane)

This document explains how to generate MOM6-readable wind/wave forcing input files for:

SCM (1D) cases forced by LES wind stress + Stokes-drift bands

Idealized 3D hurricane cases and their matched 1D sampling cases (using WW3 output)

All scripts are provided as Jupyter notebooks (.ipynb).

1. SCM (1D) forcing from LES (.mat)
What you generate

MOM6-readable forcing files containing:

Wind stress: taux, tauy

Wave/Stokes drift bands: Usx_1, Usy_1, …, Usx_k, Usy_k

Input

LES datasets stored as MATLAB files: .mat

Output (NetCDF)

The forcing NetCDF must include (minimum):

Lat, Lon

Time (unlimited dimension)

wavenumber (for wave-band information)

wind + wave variables (taux, tauy, Usx_*, Usy_*, …)

Important (FMS interpolation)

Lat, Lon, and Time must cover the full range of the target cases; otherwise FMS interpolation may error out.
A safe choice is to set:

Lat = [-10, 10]

Lon = [-10, 10]
(or any sufficiently wide range)

Notebooks

SCM_windwaveInputfiles.ipynb
Generates SCM forcing for MLD = 32 m cases (e.g., 5 m/s and 10 m/s wind–wave forcing).

SCM_windwaveInputfiles2.ipynb
Alternative generator mainly used for MLD = 10 m SCM cases (e.g., SCM vs LES comparisons at 5 m/s and 10 m/s).

2. Idealized hurricane forcing (3D + matched 1D)

Goal: Use consistent forcing between 3D and 1D cases.

2.1 Wind forcing

3D idealized hurricane wind is prescribed from the theoretical model → no external wind forcing file required for the 3D run.

We still provide wind-input generation notebooks for workflow completeness and for 1D matched cases:

Notebooks

3DwindInput.ipynb
Generates wind input files for 3D setups (includes 5 m/s and 2 m/s cases).

1DwindInput.ipynb
Generates wind files for the 1D sampling cases.

2.2 Wave forcing (WW3 → MOM6 format)

Wave forcing requires a separate 3D WW3 run (wave-only) to obtain the wave field / Stokes drift.

Workflow

Run the 3D wave-only case: 3D_wave/

Obtain WW3 output (example):
ww3.201104_usp.nc

Convert WW3 output into MOM6-readable banded Stokes drift file (example):
StokesDriftBands_201104_3DHurricaneCases.nc

Use this file as the 3D wave forcing input for MOM6

Generate matched 1D wave forcing by sampling the 3D forcing at specified times/locations

Notebook

3D1D_waveInput.ipynb

Converts WW3 output → MOM6 banded Stokes drift forcing

Produces 3D wave forcing input

Produces 1D wave forcing input by time/location sampling

3. Notebooks overview

SCM_windwaveInputfiles.ipynb — SCM forcing from LES (.mat), MLD=32 m

SCM_windwaveInputfiles2.ipynb — SCM forcing from LES (.mat), MLD=10 m

3DwindInput.ipynb — wind input generation for 3D workflows (5 m/s, 2 m/s)

1DwindInput.ipynb — wind input generation for matched 1D cases

3D1D_waveInput.ipynb — WW3 → MOM6 wave forcing (3D + 1D sampling)

hybridpaper_diagrams.ipynb — postprocessing + paper figures/diagrams

4. File management note

All wind/wave forcing files are copied into a central folder, and symbolic links are created from each case directory to the shared forcing files.

5. Case directory list (quick reference)
1D matched cases (idealized hurricane sampling)

1D_hybrid_ST
Matched SCM cases for idealized hurricane conditions using the hybrid scheme without Langmuir turbulence (ST only). 360 sampling locations.

1D_hybrid_ST+LT
Same as 1D_hybrid_ST, but with Langmuir turbulence (LT).

3D idealized hurricane cases (5 km examples)

3D_hybrid_oneway_5km
3D idealized hurricane (5 km) using hybrid scheme; LT included via one-way coupling.

3D_hybridST_5km
Same configuration but no LT. Includes grid-build scripts for 25 km / 10 km / 5 km / 1 km / 500 m.

3D_hybrid_twoway_5km
Same as 3D_hybrid_oneway_5km, but uses two-way coupling to account for LT effects.

Wave-only run (for generating forcing)

3D_wave
WW3 wave-only case used to generate wave forcing inputs for both 3D and 1D runs.

SCM configurations

SCM_GLS_ST — SCM with second-moment closure (GLS), ST only

SCM_MY_ST — SCM with second-moment closure (MY), ST only

SCM_MY_LT — SCM with second-moment closure (MY), ST+LT

SCM_hybrid — SCM using various mixing schemes (ePBL, KPP, hybrid, shear)

Northwest Atlantic (Ivan 2004)

NWA12_Ivan
1/12° NW Atlantic Ivan (2004), hybrid kappa-shear + ePBL; includes ST and one-way LT options.

NWA12_wave_Ivan
Same as NWA12_Ivan, but with two-way coupling.
