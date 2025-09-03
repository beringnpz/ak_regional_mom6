#!/bin/bash
#
# Syntax: ./calculate_surveyregion_averages.sh simname datafol
#

simname=$1
datafol=$2

#-------------------
# Averaging function
#-------------------

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

#-------------------
# Apply to files
#-------------------

# Static file with masking variables

maskfile=${simfol}/Level1-2/${simname}_ocean_static_ak.nc

# Apply to daily output files

dailyfiles=( ${simfol}/Level1-2/${simname}_selected_daily_*.nc )

for nepdailyfile in "${dailyfiles[@]}"; do
  
  svyavgfile=${nepdailyfile/Level1-2/Level3}
  svyavgfile=${svyavgfile/selected_daily/surveyregionavg}
  
  if ! test -f $svyavgfile; then
  
    echo "Calculating regional averages: ${svyavgfile}" 
   
    # Regional averages for all variables in file

    regionavg $nepdailyfile $svyavgfile $maskfile

  fi
done

# Apply to forecast files
# TODO: revisit after updates to calculate_persis_forecast

svyavgfilefc=${simfol}/Level3/${simname}_surveyregionavg_forecast_${fcyear}.nc

if ! test -f $svyavgfilefc; then
  echo "Calculating regional averages for forecast: ${fcyear}" 
   
  # Regional averages for all variables in file

  regionavg $fcfile $svyavgfilefc $maskfile
    
fi