#!/bin/bash
y_km=25_0.5h
rm -rf RESTART/*     2>/dev/null || true
rm -r 2*prog.nc ocean.stats* 2*surffluxes.nc 2*visc.nc Vertical_coordinate.nc ave_prog.nc
srun -n 100 ../../../build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6   \
    &> MOM6_${y_km}.log &   # redirect stdout/stderr if you like
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
cp -rf 20110401.prog.nc ocean.stats* 20110401.surffluxes.nc 20110401.visc.nc results/
cp -rf diag_table MOM_parameter_doc.all MOM_override MOM_input results/
cp -rf g*nml results/
mv results "./${y_km}km_results"
