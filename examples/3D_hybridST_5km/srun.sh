#!/bin/bash

module load gcp
echo "Model started:  " `date`
y_km=5_ustar0p01

rm -rf RESTART/*     2>/dev/null || true
rm -r 20*nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
current_folder=$(basename "$PWD")
JOB_ID=$(sbatch --parsable job_script.sh 1)  # First task with argument "1"
echo "Submitted Job ID: $JOB_ID"

# Monitor the job until it finishes
while true; do
    # Check if the job is still in the queue
    squeue --job $JOB_ID &> /dev/null
    JOB_EXISTS=$?  # Capture the exit status of squeue
    if [[ $JOB_EXISTS -ne 0 ]]; then
        echo " Job is no longer in the queue, breaks"
        break
    fi
    sleep 5  # Check every 5 seconds
done

# Check the job completion status using sacct
JOB_STATE=$(sacct -j $JOB_ID --format=State --noheader | tail -n 1 | awk '{print $1}')
echo "Job ID $JOB_ID has completed with status: $JOB_STATE"

# Handle job success or failure
if [[ "$JOB_STATE" == "COMPLETED" ]]; then
    echo "Job completed successfully. Proceeding..."
else
    echo "Job failed or was cancelled. Exiting..."
    exit 1
fi
mkdir results
cp -rf 20110401.prog.nc prog.nc ocean.stats* 20110401.surffluxes.nc surffluxes.nc 20110401.visc.nc visc.nc results/
cp -rf diag_table MOM_parameter_doc.all MOM_override MOM_input results/
cp -rf g*nml results/
cp -rf RESTART results/
mv results "./${y_km}_results"
#tar -czf ${y_km}_results.tar.gz ${y_km}_results
#gcp -r ${y_km}_results.tar.gz gfdl:/archive/Qian.Xiao/Qian.Xiao/MOM_GLS_LT/3D_hurricane/$current_folder/
gcp -r ${y_km}_results gfdl:/archive/Qian.Xiao/Qian.Xiao/FMS_Wave_Coupling_GOTM_kapp/3D_hurricane/$current_folder/
if [[ $? -eq 0 ]]; then
  echo "File transfer was successful."
else
  echo "File transfer failed.try to create a folder"
  mkdir $current_folder
  cp ${y_km}_results $current_folder/
  gcp -r $current_folder gfdl:/archive/Qian.Xiao/Qian.Xiao/MOM_GLS_LT/3D_hurricane/
          if [[ $? -eq 0 ]]; then
                  echo "file transfer succedded"
          else
                  echo "file transfer failed even with new filename created"
         fi
   rm -rf $current_folder
fi


