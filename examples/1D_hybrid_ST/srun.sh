#!/bin/bash

module load gcp
echo "Model started:  " `date`
start_y=2011
start_m=04
start_d=01

NUM_TOT=9  ###overall loop number=y location
rm -rf RESTART/*     2>/dev/null || true
rm -r 20*nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
rm restar*ww3 out*ww3
current_folder=$(basename "$PWD")

y_loc=(m238 m213 m188 m163 m138 m113 m88 m63 m38 m13 p13 p38 p63 p88 p113 p138 p163 p188 p213 p238)

file_name=('Wind_5mps_2998.75km_m237.5km.nc' 'Wind_5mps_2998.75km_m212.5km.nc' 'Wind_5mps_2998.75km_m187.5km.nc' 'Wind_5mps_2998.75km_m162.5km.nc' 'Wind_5mps_2998.75km_m137.5km.nc' 'Wind_5mps_2998.75km_m112.5km.nc' 'Wind_5mps_2998.75km_m87.5km.nc' 'Wind_5mps_2998.75km_m62.5km.nc' 'Wind_5mps_2998.75km_m37.5km.nc' 'Wind_5mps_2998.75km_m12.5km.nc' 'Wind_5mps_2998.75km_p12.5km.nc' 'Wind_5mps_2998.75km_p37.5km.nc' 'Wind_5mps_2998.75km_p62.5km.nc' 'Wind_5mps_2998.75km_p87.5km.nc' 'Wind_5mps_2998.75km_p112.5km.nc' 'Wind_5mps_2998.75km_p137.5km.nc' 'Wind_5mps_2998.75km_p162.5km.nc' 'Wind_5mps_2998.75km_p187.5km.nc' 'Wind_5mps_2998.75km_p212.5km.nc' 'Wind_5mps_2998.75km_p237.5km.nc')
NUM_TOT=${#y_loc[@]}

for (( i_loop=0; i_loop<NUM_TOT; i_loop++ )); do
    rm -rf RESTART/*     2>/dev/null || true
    rm restar*.ww3 out_grd.ww3
    rm -r 20110401.*.nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
  fn="${file_name[i_loop]}"
  (
    cd INPUT || { echo "ERROR: missing INPUT dir"; exit 1; }
    rm -f Wind.nc
    ln -sf "../../../forcing_1D_5mps/${fn}" Wind.nc
  )
    # get the y‐value for this iteration (will be empty if i >= length)
    y_km="${y_loc[i_loop]}"
    y_m=$(( y_km * 1000 ))
    y=$(printf "%.1E" "$y_m")
    if [[ -z "$y" ]]; then
      echo "WARNING: y_loc[$i_loop] is empty, skipping"
      continue
    fi

    echo "  ---- Iteration $i_loop: setting IDL_HURR_SCM_LOCY = $y ----"

    ../../../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6   \
        &> MOM6_${y_fmt}.log &   # redirect stdout/stderr if you like
    PID=$!
    echo "  → Started PID $PID for y=$y_fmt at $(date)"

    #  ▶ poll every 5s until it exits
    while kill -0 "$PID" 2>/dev/null; do
      echo "  → PID $PID still running… $(date)"
      sleep 5
    done

    #  ▶ now it’s gone—collect its exit code
    wait "$PID"
    ret=$?
    if [[ $ret -eq 0 ]]; then
      echo "Local run PID $PID finished SUCCESS at $(date)"
    else
      echo "Local run PID $PID finished FAILURE (exit $ret) at $(date)"
      exit 1
    fi
    mkdir results
    mv 20110401.ave_prog.nc ave_prog.nc
    mv 20110401.prog.nc prog.nc
    mv 20110401.surffluxes.nc surffluxes.nc
    mv 20110401.visc.nc visc.nc
    cp -rf prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc results/
    cp -rf MOM_IC.nc diag_table MOM_parameter_doc.all MOM_override MOM_input results/
#    cd WW3/PostProc
#        ../../../../../../build/ncrc6.intel23/ww3_ounf/REPRO/ww3_ounf > ../../log.ww3_ounf
#        if [[ -e "ww3.$start_y$start_m.nc" && -e "ww3.$start_y${start_m}_usp.nc" ]]; then
#            echo "Task 0 ww3 postprocessing profile successfully. Continuing..."
#        else
#            echo "Error in Task 0. No files ww3.$start_y$start_m.nc or ww3.$start_y${start_m}_usp.nc have been found."
#        fi
#        mv ww3.$start_y$start_m.nc ww3.$start_y$start_m$start_d.nc
#        mv ww3.${start_y}${start_m}_usp.nc ww3.${start_y}${start_m}${start_d}_usp.nc
#        cp -rf ww3*nc ../../results/
#        cp -rf ww3_ounf.inp ../../results/
#        rm ww3.$start_y*nc
#    cd ../../
    mv results "./${y_km}_results"
    #tar -czf ${y_km}_results.tar.gz ${y_km}_results
    #gcp -r ${y_km}_results.tar.gz gfdl:/archive/Qian.Xiao/Qian.Xiao/MOM_GLS_LT/3D_hurricane/$current_folder/
#    gcp -r ${y_km}_results gfdl:/archive/Qian.Xiao/Qian.Xiao/FMS_Wave_Coupling_GOTM_kapp/3D_hurricane/1DVS3D/$current_folder/
#    if [[ $? -eq 0 ]]; then
#      echo "File transfer was successful."
#    else
#      echo "File transfer failed.try to create a folder"
#      mkdir $current_folder
#      cp ${y_km}_results $current_folder/
#      gcp -r $current_folder gfdl:/archive/Qian.Xiao/Qian.Xiao/FMS_Wave_Coupling_GOTM_kapp/3D_hurricane/1DVS3D/
#              if [[ $? -eq 0 ]]; then
#                      echo "file transfer succedded"
#              else
#                      echo "file transfer failed even with new filename created"
#             fi
#       rm -rf $current_folder
#    fi

done
mkdir noLT_100
mv *_results noLT_100/
