#!/bin/bash
#
# Syntax: ./persis_forecast.sh fcyear fcdoy
#
# Where fcyear is he last partially-completed year that will be used for the persistence forecast

fcyear=$1 
fcdoy=$2
simname=$3
datafol=$4

simfol=${datafol}/${simname}

# Calculate daily forecast using persistence anomaly from last time step

anomfilefc=${simfol}/Level3/${simname}_daily_anomaly_${fcyear}.nc
fcfile=${simfol}/Level3/${simname}_forecast_${fcyear}.nc

if ! test -f $fcfile; then
  echo "Calculating persistence forecast for ${fcyear}"

  if (( $fcdoy -eq 0 )); then
    ntime1=$(cdo -ntime ${anomfilefc})
  else
    ntime1=$fcdoy
  fi
  ntime2=$(cdo -ntime ${climfile})  

  cdo -select,timestep=${ntime1} $anomfilefc tmp1.nc
  cdo -select,timestep=${ntime1}/${ntime2} $climfile tmp2.nc

  cdo -add tmp2.nc tmp1.nc $fcfile # order important, inherits year from first 
  
  rm tmp1.nc tmp2.nc
fi

# Regional averages for forecast file

svyavgfilefc=${simfol}/Level3/${simname}_surveyregionavg_forecast_${fcyear}.nc

if ! test -f $svyavgfilefc; then
  echo "Calculating regional averages for forecast: ${fcyear}" 
   
  # Regional averages for all variables in file

  regionavg $fcfile $svyavgfilefc $maskfile
    
fi