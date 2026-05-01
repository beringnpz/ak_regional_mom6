#!/bin/bash
#
# Syntax: ./calculate_clim_anom.sh fcyear
#
# Where fcyear is he last partially-completed year that will be used for the persistence forecast
# TODO: read simfol from ../simulation_data/data_folder.txt


if [ "$#" -ne 1 ]; then
    echo "Script expects one input (forecast year), exiting"
    exit 9
fi

#-----------------
# Setup
#-----------------

fcyear=$1 # Renaming for ease of reading
simname=mom6nep_hc202411
yr1=1993

simfol=../${simname}

#-----------------
# Functions
#-----------------

regionavg () {
  # $1 = input file to be averaged
  # $2 = output file
  # $3 = masking file
  
  for reg in {1..5}; do

    echo "  Region ${reg}"
  
    # Add mask variable to file
  
    cp $1 tmpo.nc
    ncks -A -v areacello,mask_survey $3 tmpo.nc
  
    # Add cold pool variables

    cdo -expr,'cpool2p0=tob<2.0; cpool1p0=tob<1.0; cpool0p0=tob<0.0' $1 tmpcp.nc
    ncks -A tmpcp.nc tmpo.nc

    # Calculate weighted average
  
    ncwa -w areacello -B "mask_survey = ${reg}" -a xh,yh tmpo.nc tmpreg${reg}.nc
  
    # Clean up 
  
    rm tmpo.nc tmpcp.nc
  done
  
  # Concatenate regions along new dimension
  
  ncecat -u surveyregion tmpreg*.nc $2
  rm tmpreg*.nc

}

#-----------------
# Calculations
#-----------------

# # Add omega variable

# for (( yr=$yr1; yr<=$fcyear; yr++ )); do

#   nepdailyfile=${simfol}/Level1-2/${simname}_selected_daily_${yr}.nc

#   if (cdo -showname $nepdailyfile | grep -q btm_omega_arag); then

#     cdo -expr,'btm_omega_arag=btm_co3_ion/btm_co3_sol_arag' ${simfol}/Level1-2/${simname}_selected_daily_${yr}.nc tmp.nc

#   fi
# done


# Calculate daily climatology for 30-year period (1993-2022) 
# (TODO: Need to figure out how to dynamically build that glob in a way cdo accepts)

climfile=${simfol}/Level3/${simname}_daily_clim_1993-2022.nc

if ! test -f $climfile; then
  echo "Building 1993-2022 climatology"
  cdo -ydaymean -selyear,1993/2022 -cat '../mom6nep_hc202411/Level1-2/mom6nep_hc202411_selected_daily*.nc' $climfile
  # cdo -ydaymean -selyear,1993/2022 -cat "${simfol}/Level1-2/$simname_selected_daily*.nc" $climfile
fi

# Calculate daily anomaly relative to climatology

for (( yr=$yr1; yr<=$fcyear; yr++ )); do
  
  nepdailyfile=${simfol}/Level1-2/${simname}_selected_daily_${yr}.nc
  anomfile=${simfol}/Level3/${simname}_daily_anomaly_${yr}.nc

  if ! test -f $anomfile; then
    echo "Calculating anomalies: ${yr}"
    cdo ydaysub $nepdailyfile $climfile $anomfile
  fi

done

# Calculate daily forecast using persistence anomaly from last time step

anomfilefc=${simfol}/Level3/${simname}_daily_anomaly_${fcyear}.nc
fcfile=${simfol}/Level3/${simname}_forecast_${fcyear}.nc

if ! test -f $fcfile; then
  echo "Calculating persistence forecast for ${fcyear}"

  ntime1=$(cdo -ntime ${anomfilefc})
  ntime2=$(cdo -ntime ${climfile})  

  cdo -select,timestep=${ntime1} $anomfilefc tmp1.nc
  cdo -select,timestep=${ntime1}/${ntime2} $climfile tmp2.nc

  cdo -add tmp2.nc tmp1.nc $fcfile # order important, inherits year from first 
  
  rm tmp1.nc tmp2.nc
fi

# Regionally-averaged timeseries based on AK survey regions

for (( yr=$yr1; yr<=$fcyear; yr++ )); do
  
  nepdailyfile=${simfol}/Level1-2/${simname}_selected_daily_${yr}.nc
  svyavgfile=${simfol}/Level3/${simname}_surveyregionavg_${yr}.nc
  maskfile=${simfol}/Level1-2/ocean_static_ak.nc
  
  if ! test -f $svyavgfile; then
  
    echo "Calculating regional averages: ${yr}" 
   
    # Regional averages for all variables in file

    regionavg $nepdailyfile $svyavgfile $maskfile

  fi
done

# Regional averages for forecast file

svyavgfilefc=${simfol}/Level3/${simname}_surveyregionavg_forecast_${fcyear}.nc

if ! test -f $svyavgfilefc; then
  echo "Calculating regional averages for forecast: ${fcyear}" 
   
  # Regional averages for all variables in file

  regionavg $fcfile $svyavgfilefc $maskfile
    
fi


