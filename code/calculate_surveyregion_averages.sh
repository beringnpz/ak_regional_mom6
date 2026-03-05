#!/bin/bash
#
# Syntax: ./calculate_surveyregion_averages.sh simname datafol
#

# simname=$1
# datafol=$2

#-----------------
# Setup
#-----------------

source cefi_name_tables.sh

ppfol="."
exptype="hcast"

# Input parsing

while [[ $# -gt 0 ]]; do
  case $1 in
    --ppdir)
      ppfol=$2
      shift # past argument
      shift # past value
      ;;
    --region)
      region=$2
      shift
      shift
      ;;
    --subdomain)
      subdomainstr=$2
      shift
      shift
      ;;
    --exptype)
      exptype=$2
      shift
      shift
      ;;
    --release)
      release=$2
      shift
      shift
      ;;
    --maskbase)
      maskbase=$2
      shift
      shift
      ;;
    --maskdatestr)
      maskdatestr=$2
      shift
      shift
      ;;
    -h|--help)
       echo "$USEAGE"
       exit
       ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
  esac
done

if [[ "${subdomainstr}" == "iq*" ]] ; then
    subdomainstr_long=${tbl_subdomain[${subdomainstr}]}
else
    subdomainstr_long=${subdomainstr}
fi


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

# Base folder

basefol="${ppfol}/${tbl_region[${region}]}/${subdomainstr_long}/${tbl_exptype[${exptype}]}/"
# daily/extra/${release}"

# Static file with masking variables

maskfile=${basefol}/static/extra/${maskbase}.${region}.${subdomainstr}.${exptype}.static.${release}.${maskdatestr}.nc

# Apply to daily output, forecast, anomaly, and climatology files

files = (${basefol}/daily/raw/*.nc
         ${basefol}/daily/extra/fcpersis_*.nc
         ${basefol}/daily/extra/anom_*.nc
         ${basefol}/daily/extra/clim*.nc
        )

# Output folder (considering survey region average indices to be its own 
# subdomain, since this is somewhat independent of exact subdomain)

outfol="${ppfol}/${tbl_region[${region}]}/ak_surveyregion_avg/${tbl_exptype[${exptype}]}/daily/extra"

# files = (${simfol}/Level1-2/${simname}_selected_daily_*.nc
#          ${simfol}/Level3/${simname}_forecast_*.nc
#          ${simfol}/Level3/${simname}_daily_anomaly*.nc
#          ${simfol}/Level3/${simname}_daily_clim*.nc
#         )

# dailyfiles=( ${simfol}/Level1-2/${simname}_selected_daily_*.nc )

if [ ! -d "${outfol}" ]; then
    mkdir -p "${outfol}"
fi

for fname in "${files[@]}"; do

  svyavgfile=$( basename "${fname}" )
  svgavgfile=${outfol}/${svgavgfile/${subdomainstr}/aksvyreg}

  # svyavgfile=${fname/Level1-2/Level3}
  # svgavgfile=${svyavgfile/Level3/Level3\/surveyregionavg}
  # svyavgfile=${svyavgfile/.nc/.svyreg.nc}
  
  if ! test -f $svyavgfile; then
  
    echo "Calculating regional averages: ${svyavgfile}" 
   
    # Regional averages for all variables in file

    regionavg $fname $svyavgfile $maskfile

  fi
done

# mv ${simfol}/Level3/*.svyreg.nc 

# # Apply to forecast files
# # TODO: revisit after updates to calculate_persis_forecast

# fcfiles=( ${simfol}/Level3/${simname}_forecast_*.nc )

# # svyavgfilefc=${simfol}/Level3/${simname}__*.nc)

# if ! test -f $svyavgfilefc; then
#   echo "Calculating regional averages for forecast: ${fcyear}" 
   
#   # Regional averages for all variables in file

#   regionavg $fcfile $svyavgfilefc $maskfile
    
# fi