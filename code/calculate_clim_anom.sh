#!/bin/bash
#
# Syntax: ./calculate_clim_anom.sh simname datafol
#


if [ "$#" -ne 2 ]; then
    echo "Script expects 2 inputs (simname, datafol), exiting"
    exit 9
fi

#-----------------
# Setup
#-----------------

simname=$1
datafol=$2

yr1=1993

simfol=${datafol}/${simname}

#-----------------
# Calculations
#-----------------

# Calculate daily climatology for 30-year period (1993-2022) 
# (TODO: Need to figure out how to dynamically build that glob in a way cdo accepts)

climfile=${simfol}/Level3/${simname}_daily_clim_1993-2022.nc

#globstr="\'${simfol}/Level1-2/${simname}_selected_daily*.nc\'"

if ! test -f $climfile; then
  echo "Building 1993-2022 climatology"
  cdo -ydaymean -selyear,1993/2022 -cat '/gscratch/cicoes/GR011846_reem/CEFI_data/mom6nep_hc202507/Level1-2/mom6nep_hc202507_selected_daily*.nc' $climfile
  #cdo -ydaymean -selyear,1993/2022 -cat $globstr $climfile
fi

# Calculate daily anomaly relative to climatology

dailyfiles=( ${simfol}/Level1-2/${simname}_selected_daily_*.nc )

for nepdailyfile in "${dailyfiles[@]}"; do
  
  anomfile=${nepdailyfile/Level1-2/Level3}
  anomfile=${anomfile/selected_daily/daily_anomaly}

  if ! test -f $anomfile; then
    echo "Calculating anomalies: ${anomfile}"
    cdo ydaysub $nepdailyfile $climfile $anomfile
  fi

done

