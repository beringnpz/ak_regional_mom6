#!/bin/bash
#
# Syntax: ./calculate_clim_anom.sh simname datafol
#
# Where fcyear is he last partially-completed year that will be used for the persistence forecast
# TODO: read simfol from ../simulation_data/data_folder.txt


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
# Functions
#-----------------

regionavg () {
  # $1 = input file to be averaged
  # $2 = output file
  # $3 = masking file
  
  regnum=(6 47 52 78 98 143)

  for reg in "${regnum[@]}"; do

    echo "  Region ${reg}"
  
    # Add mask variable to file

    cp $1 tmpo.nc
    ncks -A -v areacello,mask_survey_area $3 tmpo.nc
  
    # Add cold pool variables
 
    cdo -expr,'cpool2p0=tob<2.0; cpool1p0=tob<1.0; cpool0p0=tob<0.0' $1 tmpcp.nc
    ncks -A -v areacello,mask_survey_area $3 tmpcp.nc

    # Calculate weighted average

    ncwa    -w areacello -B "mask_survey_area = ${reg}" -a ih,jh tmpo.nc  tmpreg${reg}.nc
    ncwa -A -w areacello -B "mask_survey_area = ${reg}" -a ih,jh tmpcp.nc tmpreg${reg}.nc
  
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

# Regionally-averaged timeseries based on AK survey regions

maskfile=${simfol}/Level1-2/${simname}_ocean_static_ak.nc

for nepdailyfile in "${dailyfiles[@]}"; do
  
  svyavgfile=${nepdailyfile/Level1-2/Level3}
  svyavgfile=${svyavgfile/selected_daily/surveyregionavg}
  
  if ! test -f $svyavgfile; then
  
    echo "Calculating regional averages: ${svyavgfile}" 
   
    # Regional averages for all variables in file

    regionavg $nepdailyfile $svyavgfile $maskfile

  fi
done

