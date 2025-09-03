#!/bin/bash
#
# Syntax: ./persis_forecast.sh filedatestr tindex simname datafol
#
# Where fcyear is he last partially-completed year that will be used for the persistence forecast

if [ "$#" -ne 4 ]; then
    echo "Script expects 4 inputs (filedatestr, tindex, simname, datafol), exiting"
    exit 9
fi

#-----------------
# Setup
#-----------------

filedatestr=$1 
tindex=$2
simname=$3
datafol=$4

simfol=${datafol}/${simname}

#-----------------
# Calculations
#-----------------

# Determine time step to use, and format date strings accordingly

anomfilefc=${simfol}/Level3/${simname}_daily_anomaly_${filedatestr}.nc

if (( $tindex -eq 0 )); then
    ntime1=$(cdo -ntime ${anomfilefc})
  else
    ntime1=$fcdoy
fi

fcdatestr=$(cdo -showdate -select,timestep=${ntime1} $anomfilefc)
fcdatestr=${fcdatestr//-/}

fcfile=${simfol}/Level3/${simname}_forecast_${fcdatestr}.nc

# Calculate daily forecast using persistence anomaly from one time step

climfile=${simfol}/Level3/${simname}_daily_clim_1993-2022.nc

# Calculate daily forecast using persistence anomaly from one time step 
# added to climatology, from selected time step until the end of the 
# calendar year

if ! test -f $fcfile; then
  echo "Calculating persistence forecast for ${fcyear}"

  ntime2=$(cdo -ntime ${climfile})  

  cdo -select,timestep=${ntime1} $anomfilefc tmp1.nc
  cdo -select,timestep=${ntime1}/${ntime2} $climfile tmp2.nc

  cdo -add tmp2.nc tmp1.nc $fcfile # order important, inherits year from first 
  
  rm tmp1.nc tmp2.nc
fi