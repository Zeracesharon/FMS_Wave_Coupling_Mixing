#!/bin/bash

module load gcp
echo "Model started:  " `date`
start_y=2011
start_m=04
start_d=01

rm -rf RESTART/*     2>/dev/null || true
rm -r 20*nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
rm restar*ww3 out*ww3
current_folder=$(basename "$PWD")
Y0_KM=900.0
DYPRIME_START=-897.5
DYPRIME_STEP=5.0
NUM_TOT=360

for (( i_loop=1; i_loop<=NUM_TOT; i_loop++ )); do
    rm -rf RESTART/*     2>/dev/null || true
    rm -f restar*.ww3 out_grd.ww3
    rm -f 20110401.*.nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc

    sid=$(printf "%03d" "$i_loop")
    fn="Wind_5mps_${sid}.nc"

    dyprime_km=$(awk -v i="$i_loop" -v a="$DYPRIME_START" -v s="$DYPRIME_STEP" 'BEGIN{printf "%.1f", a + s*(i-1)}')
    y_km=$(awk -v y0="$Y0_KM" -v dy="$dyprime_km" 'BEGIN{printf "%.1f", y0 + dy}')

    tag=$(awk -v dy="$dyprime_km" 'BEGIN{
        if (dy < 0) printf "m%.1fkm", -dy;
        else        printf "p%.1fkm", dy;
    }')

    y_m=$(awk -v y="$y_km" 'BEGIN{printf "%.1f", y*1000.0}')
    y=$(awk -v ym="$y_m"  'BEGIN{printf "%.1E", ym}')

    echo "---- Iteration ${sid}: dyprime=${dyprime_km} km  y=${y_km} km  (IDL_HURR_SCM_LOCY=${y}) ----"

    # link forcing (WIND + WAVE)
    (
      cd INPUT || { echo "ERROR: missing INPUT dir"; exit 1; }

      # wind
      rm -f Wind.nc
      ln -sf "../../../forcing_1D_5mps/Wind_5mps_${sid}.nc" Wind.nc

      # wave (MOM6 override-format bands file for this station)
      rm -f StokesDriftBands_2004.nc
      ln -sf "../../../1D_wave_forcing_5mps/Wave_5mps_${sid}.nc" StokesDriftBands_2004.nc
    )

    ../../../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6 \
        &> "MOM6_${sid}_${tag}.log" &
    PID=$!

    while kill -0 "$PID" 2>/dev/null; do
      sleep 5
    done
    wait "$PID"; ret=$?
    if [[ $ret -ne 0 ]]; then
      echo "Run FAILURE (exit $ret) sid=${sid} tag=${tag}"
      exit 1
    fi

    mkdir -p results
    mv 20110401.ave_prog.nc ave_prog.nc
    mv 20110401.prog.nc prog.nc
    mv 20110401.surffluxes.nc surffluxes.nc
    mv 20110401.visc.nc visc.nc

    cp -rf prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc results/
    cp -rf MOM_IC.nc diag_table MOM_parameter_doc.all MOM_override MOM_input results/

    mv results "./${sid}_${tag}_results"
done
mkdir oneway_all_360
mv *_results oneway_all_360/
