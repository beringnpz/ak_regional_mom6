#!/bin/bash
#
# Syntax: ./calculate_clim_anom.sh targetyear
#
# Where target year is he last partially-completed year that will be used for the persistence forecast

if [ "$#" -ne 1 ]; then
    echo "Script expects one input (target year), exiting"
    exit 9
fi

#-----------------
# Functions
#-----------------

regionavg () {
  # $1 = input file to be averaged
  # $2 = output file
  # $3 = masking file
  
  for reg in {1..5}; do

    echo "  Region ${reg}"
   
    # Put ice on xh,yh grid for ease of calculations

    ncks -O -v    siconc,xT,yT $1 tmpi.nc # ice grid vars
    ncks -O -x -v siconc,xT,yT $1 tmpo.nc # ocean grid vars
  
    ncrename -d xT,xh -d yT,yh tmpi.nc
    ncks -A tmpi.nc tmpo.nc
  
    # Add mask variable to file
  
    ncks -A -v areacello,mask_survey $3 tmpo.nc
  
    # Calculate weighted average
  
    ncwa -w areacello -B "mask_survey = ${reg}" -a xh,yh tmpo.nc tmpreg${reg}.nc
  
    # Clean up 
  
    rm tmpi.nc tmpo.nc
  done
  
  # Concatenate regions along new dimension
  
  ncecat -u surveyregion tmpreg*.nc $2
  
  rm tmpreg*.nc

}

#-----------------
# Setup
#-----------------

climfile=./pp/nep202411_daily_clim_1993-2022.nc

# Calculate daily climatology for 30-year period (1993-2022) 

if ! test -f $climfile; then
  echo "Building 1993-2022 climatology"
  cdo -ydaymean -selyear,1993/2022 -cat 'mom6nep_selected/*.nc' $climfile
fi

# Calculate daily anomaly relative to climatology

yr1=1993

for (( yr=$yr1; yr<=$1; yr++ )); do
  
  nepdailyfile=./mom6nep_selected/nep202411_selected_daily_${yr}.nc
  anomfile=./pp/nep202411_daily_anomaly_${yr}.nc

  if ! test -f $anomfile; then
    echo "Calculating anomalies: ${yr}"
    cdo ydaysub $nepdailyfile $climfile $anomfile
  fi

done

# Calculate daily forecast using persistence anomaly from last time step

anomfiletarget=./pp/nep202411_daily_anomaly_${1}.nc
fcfile=./pp/nep202411_forecast_${1}.nc


if ! test -f $anomfiletarget; then
  echo "Calculating persistence forecast for ${1}"

  ntime1=$(cdo -ntime ${anomfiletarget})
  ntime2=$(cdo -ntime ${climfile})  

  cdo -select,timestep=${ntime1} $anomfiletarget tmp1.nc
  cdo -select,timestep=${ntime1}/${ntime2} $climfile tmp2.nc

  cdo -add tmp1.nc tmp2.nc $fcfile
  
  rm tmp1.nc tmp2.nc
fi

# Regionally-averaged timeseries based on AK survey regions

for (( yr=$yr1; yr<=$1; yr++ )); do
  
  nepdailyfile=./mom6nep_selected/nep202411_selected_daily_${yr}.nc
  svyavgfile=./pp/nep202411_surveyregionavg_${yr}.nc
  maskfile=./ocean_static_ak.nc
  
  if ! test -f $svyavgfile; then
  
    echo "Calculating regional averages: ${yr}" 
   
    regionavg $nepdailyfile $svyavgfile $maskfile
    
    # for reg in {1..5}; do
    #
    #   echo "  Region ${reg}"
    #
    #   # Put ice on xh,yh grid for ease of calculations
    #
    #   ncks -v    siconc,xT,yT ./mom6nep_selected/nep202411_selected_daily_${yr}.nc tmpi.nc # ice grid vars
    #   ncks -x -v siconc,xT,yT ./mom6nep_selected/nep202411_selected_daily_${yr}.nc tmpo.nc # ocean grid vars
    #
    #   ncrename -d xT,xh -d yT,yh tmpi.nc
    #   ncks -A tmpi.nc tmpo.nc
    #
    #   # Add mask variable to file
    #
    #   ncks -A -v areacello,mask_survey ocean_static_ak.nc tmpo.nc
    #
    #   # Calculate weighted average
    #
    #   ncwa -w areacello -B "mask_survey = ${reg}" -a xh,yh tmpo.nc tmp${yr}reg${reg}.nc
    #
    #   # Clean up
    #
    #   rm tmpi.nc tmpo.nc
    # done
    #
    # # Concatenate regions along new dimension
    #
    # ncecat -u surveyregion tmp${yr}reg*.nc pp/nep202411_surveyregionavg_${yr}.nc
    #
    # rm tmp${yr}reg*.nc
    
  fi
done





#     # Regional average per file
#     cp ./mom6nep_selected/nep202411_selected_daily_${yr}.nc tmp.nc
#     # Put ice on xh,yh grid for ease of calcs
#     ncks -v siconc,xT,yT tmp.nc tmpi.nc
#     ncks -O -x -v siconc,xT,yT tmp.nc tmp.nc
#     ncrename -d xT,xh -d yT,yh tmpi.nc
#     ncks -A tmpi.nc tmp.nc
#     # Add mask variable to file
#     ncks -A -v areacello,mask_survey ocean_static_ak.nc tmp.nc
#     ncwa -w areacello -B "mask_survey = ${reg}" -a xh,yh tmp.nc tmp${yr}reg${reg}.nc
#     rm tmp.nc tmpi.nc
#   done
#   ncrcat tmp*reg${reg}.nc reg${reg}.nc
#   rm tmp*reg${reg}.nc
# done
# # Concatenate regions
# ncecat -u surveyregion reg*.nc pp/nep202411_surveyregionavg.nc
# rm reg*.nc
  


