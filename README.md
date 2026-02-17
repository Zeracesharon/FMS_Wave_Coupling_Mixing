Overview
This repo adds two more vertical mixing schemes on top of the repo FMS_Wave_Coupling_Dec2020: Second momentum closure (SMC) from GOTM and a Hybrid scheme that couples kappa-shear and ePBL. Both of them can be used to account for Langmuir turbulence (LT). The installation method has been documented in FMS_Wave_Coupling_Dec2020 and will not be repeated here. Note that the current installation environment files can be used for both C5 and C6 systems in Gaea. The MOM6 source file package in this wave coupling repo is not the most updated MOM6-example so the installation will be slightly different from the most updated MOM6 guidance. This repo gives everything needed for the paper (Xiao & Reichl, 2026)-- which will be online after publication. Hereafter, the paper refers to Xiao & Reichl (2026).
 
Obtaining Wave Information (Three Options)
Generally, there are three different ways to obtain wave information: theoretical waves, wave-forcing input file from pre-run wave models (one-way coupling) and from a coupled wave model (two-way coupling).

1) Theoretical Waves (LF17)
   
1.	The theoretical wave aligns with the wind direction and uses wind information to directly calculate wave properties (wavenumber and wave height), and uses exponential decay for the vertical profile, and a finite difference method to obtain the Stokes drift shear. To use theoretical waves for LT, set USE_WAVES = True; WAVE_METHOD = LF17.
   
2) One-way Coupling (Wave Forcing Input File)
1.	For one way coupling, the wave information comes from a wave input file that contains dimensions of Time (unlimited); Lat and Lon; and wave number (we typically use 3-14 different wavenumber components ranging from 0 to 4). Variables are time, lon, lat and Usx_k (the x-direction of surface Stokes drift for a given wavenumber) and Usy_k (the y-direction of surface Stokes drift for a given wavenumber) where the time step and spatial domain should cover the required region and time duration otherwise FMS interpolation will report error. The recommended value and number of wavenumbers will be denoted in specific examples.
   
Recommended wavenumbers used in the paper
We use 14 wavenumbers (WAVENUMBERS = [0.006, 0.01 , 0.02 , 0.04 , 0.06 , 0.08 , 0.1 , 0.2 , 0.4 , 0.6, 0.8 , 1. , 2. , 4. ]) in our one-way coupling simulation cases in the paper. Then the vertical profile of the Stokes drift value will be obtained by assuming exponential decay for each component. The way and script to obtain the wave-forcing input file will be given in the example 3D hurricane cases.
Data table configuration (FMS)

When the wave-input files are ready, the data_table should be set as:
"OCN" , "Usx1" , "Usx1" , "INPUT/Hurr.nc" , "bilinear" , 1.0
"OCN" , "Usy1" , "Usy1" , "INPUT/Hurr.nc" , "bilinear" , 1.0
…
"OCN" , "Usx14" , "Usx14" , "INPUT/Hurr.nc" , "bilinear" , 1.0
"OCN" , "Usy14" , "Usy14" , "INPUT/Hurr.nc" , "bilinear" , 1.0
Where "INPUT/Hurr.nc" indicates the wave-forcing file.
MOM_override configuration
In MOM_override, set: USE_WAVES = True; WAVE_METHOD = SURFACE_BANDS; SURFBAND_SOURCE = DATAOVERRIDE;
SURFBAND_FILENAME = INPUT/Hurr.nc.

3) Two-way Coupling (Coupled Ocean–Wave Model)
   
2.	From a coupled model for two-way coupling, see the original folder (FMS_Wave_Coupling_Dec2020) for wave-current information exchange.
Settings
For two-way coupling setup: In input.nml set do_waves=.true., set up WW3 grid and ww3_multi.inp, in MOM_override set: USE_WAVES = True; WAVE_METHOD = SURFACE_BANDS; SURFBAND_SOURCE = COUPLER; STK_BAND_COUPLER = 3; The details of how to set up the wave model can be referred to the specific example.
If you consider the LT effects within the ePBL framework, then turn wave off, you set:
#override USE_LA_LI2016 = False; #override EPBL_LT = False; #override USE_WAVES = False
Turn waves on, you have:
#override USE_WAVES = True; #override USE_LA_LI2016 = True; #override EPBL_LT = False for theoretical waves,
#override USE_WAVES = True; #override USE_LA_LI2016 = False; #override EPBL_LT = True, for real waves (one-way+two-way)
 
Added Vertical Mixing Approaches: Overall Illustration
Here I will give an overall illustration of the two added vertical mixing approach setting and will refer to the example folder for the detailed setting. The main purpose of this is to illustrate how to set up the input files and work flow to reproduce the cases in the paper.
 
GOTM Coupling and Case Setup
1) Using SMC in GOTM
•	Use SMC in GOTM as the vertical mixing scheme. There are two relevant papers that introduce the SMC that I use in the examples:
"A generic length-scale equation for geophysical turbulence models" and "A second-moment closure model of Langmuir turbulence",
where the first one introduces GLS and the second one introduces how to use GOTM scheme choices to account for Langmuir turbulence.
Example folders
The GLS example is SCM_GLS_ST (no LT included) and the one with LT case setting is SCM_MY_LT.
Turning off overlapping schemes (to avoid redundancy)
When GOTM package is used, ePBL, KPP and kappa-shear should be closed to avoid redundancy. To set this, in MOM_override, set: USE_JACKSON_PARAM = False; ENERGETICS_SFC_PBL = False, USE_KPP = False; USE_GOTM = True, GOTM_VERTEX_SHEAR = False (Or True both works ).
What input.nml vs gotmturb.nml controls
For each case, the input.nml only controls FMS coupler with Ocean and wave, the GOTM package is controlled by gotmturb.nml, where the detailed documentation of each namelist can be refered to GOTM guidance: https://gotm.net/portfolio/ and its source code.
GOTM turbulence model choices used here (examples)
Here we only give some combination of the method that I have been used: for &turbulence: 3,2,10,3 uses GLS model for the second-order model, we chose weak equilibrium scnd=2. There are several parameters that you can play with when you want to tune the model in &turb_param: compute_kappa (True or false), kappa number (0.4), compute_c3 (true or false), Ri_st (0.25-0.38), length_lim (true or false) and galp.
&turbulence: 3,3,9,1 is the dynamic Mellor-Yamada q^2l-equation, correspondingly, you can tune the parameter in namelist &my: sq, sl, new_constr.
Adding LT on top of tuned turbulence settings

For examples that include LT, this should be on top of the model that has been tuned above. To consider LT effects, you will need wave information and modify the MOM_override to let MOM know the wave information which can be prepared and set as indicated above.
Then in GLS, the LT related terms should be included in the namelist: &generic: cpsix= 1.44, cpsi4= 0.3,
&keps: cex=1.44;ce4=0.0; &kw: cwx = 0.555, cw4 = 0.15.
In Harcourt (2013), they recommend using MY model (3,3,9,1), then the related term is &my: e6 = 7.0. The namelist &scnd: scnd_method = 4 should be changed for every turbulent model chosen in order to consider LT effects. There are some instructions in MOM_override files in each example folder too.
 
Hybrid Mixing Scheme (kappa-shear + ePBL)
2) Detailed explanation of settings
•	For the hybrid mixing scheme, the idea has been illustrated in the paper. Here we give the detailed explanation of each setting:
When to use this hybrid scheme
The hybrid mixing is designed for hurricane conditions or wind-wave forcing that includes a hurricane. If you are facing normal low to moderate wind, you can simply set EPBL_USE_SHEAR_ENERGY_COUPLING = False. Then the code will return back to the maximum combination of ePBL+kappa-shear.
Enabling hybrid mixing (required base switches)
If you want hybrid mixing, set this to be true. Note that in order to ensure this hybrid method works, ePBL and kappa-shear should all be turned on: USE_JACKSON_PARAM = True; ENERGETICS_SFC_PBL = True, USE_KPP = False; USE_GOTM = False. The hybrid energy scheme will only work when both ePBL and kappa-shear are turned on, otherwise, it will do nothing regardless of the following settings.
Option descriptions (as used in this repo)

EPBL_SHEAR_ENERGY_SHUTDOWN_ML: this decides how kappa-shear and energy-augmented ePBL combines for turbulent coefficient. If true, means above the mixed layer depth (which is identified either by using iterated MLD from ePBL or the depth of the maximum vertical temperature gradient), the energy-augmented ePBL will be the only one that controls the turbulent coefficient. If false, the maximum combination of energy-augmented ePBL and kappa-shear remains. According to my experience, this does not change the results a lot as the energy-augmented ePBL produces similar results with kappa-shear within the MLD layer.

EPBL_SHEAR_ENERGY_WARN_IF_NO_KS: this can always set to be true to warn if kappa-shear has not been used but this energy coupling hybrid method has been turned on.

EPBL_LT_RESCALE_BY_SHEAR: if you want rescaled LT effects, this needs to be true, otherwise it uses the mLT directly from RL19 (Reichl & Li, 2019).
SHEARENERGY_HBL_EST: this indicates the mixed layer depth we use to calculate the mechanical energy using shear-diagnosed buoyancy flux. True if the depth of maximum vertical temperature gradient is used, otherwise, MLD diagnosed by ePBL through energy consistent iteration is used.
USTAR_MIN_THRE: the minimum ustar value to initialize the shear-coupling. Any ustar value smaller than this, the shear-coupling workflow will be shutdown.

LT_SHEAR_FRAC_CAP; SHEAR_EPBL_CAP, SHEAR_EPBL_LT_CAP: these three cap are used to bound the ratio of enhancement on mechanical energy corresponding to shear turbulent and LT. SHEAR_EPBL_CAP confines the ratio of shear-diagnosed mechanical energy used by ePBL to its parameterized value obtained by RH18 (Reichl & Hallberg, 2018). SHEAR_EPBL_LT_CAP confines the ratio of rescaled m_LT to the value by RL19. LT_SHEAR_FRAC_CAP indicates another way to confine the rescaling factor (make sure the rescaled m_LT should be smaller than the cap multiply by the m of the baseline after augmentation). This is normally set to be very large to disable this confinement to allow SHEAR_EPBL_LT_CAP to be effective.

LT_ENHANCE_TOP is a debugging choice to allow the user to set the rescaling factor manually. Note that this works on top of the rescaling factor calculated based on rescaling schemes. To avoid double rescaling, either set SHEAR_EPBL_LT_CAP=1.0 or LT_ENHANCE_TOP=1.0. These two LT rescaling ways are effective only when EPBL_LT_RESCALE_BY_SHEAR is true.

SHEAR_LT_MODE: we chose preserve_fraction as the quantity to consider LT enhancement and rescaling. The other option to rescale is energy_gamma which rescales m*_LT by (E_used/E_ePBL)^gamma where E_used is the shear-diagnosed mechanical energy.

EPBL_LT_RESCALE_GAMMA: it can either be the gamma indicated above in energy_gamma or in preserve_fraction works as (target_full/mLTtmp)^gamma where (target_full/mLTtmp) is the calculated rescaling factor.

MIX_LEN_EXPONENT: is the exponential parameter to decide the shape of the turbulent coefficients and thus the distribution of the energy. Normally can be any value between 1 to 3.

IS_KS_WEIGHT: indicates whether to smooth the turbulent coefficient profile when EPBL_SHEAR_ENERGY_SHUTDOWN_ML is true. As you can imagine, at the depth near the diagnosed mixed layer, above it, we use the shear-energy augmented coefficient, below it we use kappa-shear. Sometimes there will be some sharp changes. To make the profile smooth, we apply a smooth weighted function at the last W_RATIO depth, we apply the smooth function that changes exponentially with depth to achieve a smooth transition.

W_RATIO: the ratio of depth that started to apply smooth function to diagnosed mixed-layer depth.

SHEAR_ALPHA: this is a debug option, the recommended value is 1. This is the ratio of the depth used to integrate the shear-diagnosed buoyancy flux to obtain mechanical energy to the mixed-layer depth diagnosed either by ePBL iteration or maximum temperature gradient. This influences how much shear energy has been transferred to augment the ePBL entrainment energy constraint.

STOKES_VF= True, STOKES_DDT=True, STOKES_PGF=True these three options controls how many Stokes drift associated terms in the wave-averaged equation in the paper are calculated in the real momentum equation.
 
Examples and Postprocessing
For each example case folder, there will be a separate documentation providing the details of setting the case, preparing the input files and the workflow of running the model. There will also be a folder that contains postprocessing file and code for the diagrams drawn in the paper.
 
Notes on MOM6 Version and Portability
As this is a two-way coupled ocean-wave model, the MOM6 version is quite old. If you want to use GOTM or hybrid approach with one-way coupling or cases without wave information, readers are referred to the package of MOM6_Mixing for the updated MOM6 version. (The corresponding MOM6 is the version of 2025 Sep). If the reader wants the most updated MOM6 version with GOTM or hybrid method, the code incorporated with these two vertical mixing schemes is easy to transplant.
 
References

Umlauf, L. and Burchard, H., 2003. A generic length-scale equation for geophysical turbulence models.
Harcourt, R.R., 2013. A second-moment closure model of Langmuir turbulence. Journal of Physical Oceanography, 43(4), pp.673-697.
Reichl, B.G. and Li, Q., 2019. A parameterization with a constrained potential energy conversion rate of vertical mixing due to Langmuir turbulence. Journal of Physical Oceanography, 49(11), pp.2935-2959.
Reichl, B.G. and Hallberg, R., 2018. A simplified energetics based planetary boundary layer (ePBL) approach for ocean climate simulations. Ocean Modelling, 132, pp.112-129.
Xiao, Q. and Reichl, B.G., 2026. A hybrid $\kappa$-shear--ePBL vertical mixing scheme with a Langmuir turbulence parameterization for tropical cyclones in a coupled ocean-wave model.

