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
    --vars)
      varstr=$2
      shift
      shift
      ;;
    --tindex)
      tindex=$2
      shift
      shift
      ;;
    --filedatestr)
      filedatestr=$2
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


# filedatestr=$1 
# # tindex=$2
# simname=$3
# datafol=$4

# simfol=${datafol}/${simname}

#-----------------
# Calculations
#-----------------

# Base folder

basefol="${ppfol}/${tbl_region[${region}]}/${subdomainstr_long}/${tbl_exptype[${exptype}]}/daily/extra/${release}"

# Parse variables

if [[ "${varstr[$i]}" == *","*  ]]; then # if multiple (contains comma)
  IFS=',' read -ra varnames <<< "${varstr[$i]}"
else
  varnames=( "${varstr[$i]}" )
fi

# Loop over variables...

for vv in "${varnames[@]}"; do
  
  # Anomaly file to use as forecast base

  anomfilefc="${basefol}/anom_${vv}.${region}.${subdomainstr}.${exptype}.daily.${release}.${filedatestr}.nc"

  # Determine time step to use, and format date strings accordingly

  if (( $tindex == 0 )); then
      ntime1=$(cdo -ntime ${anomfilefc})
    else
      ntime1=$tindex
  fi

  # Translate time index to a file start-date string

  fcdatestr=$(cdo -showdate -select,timestep=${ntime1} $anomfilefc)
  fcdatestr=${fcdatestr//-/} # remove hyphens ...
  fcdatestr=${fcdatestr// /} # ... and spaces

  # Forecast file name

  fcfile="${basefol}/fcpersist_${vv}.${region}.${subdomainstr}.${exptype}.daily.${release}.${fcdatestr}.nc"

  # Calculate daily forecast using persistence anomaly from one time step

  climfile="${basefol}/clim_${vv}.${region}.${subdomainstr}.${exptype}.daily.${release}.1993-2022.nc"

  # Calculate daily forecast using persistence anomaly from one time step 
  # added to climatology, from selected time step until the end of the 
  # calendar year

  if [ ! -f $fcfile ]; then

    ntime2=$(cdo -ntime ${climfile})  

    cdo -select,timestep=${ntime1} $anomfilefc tmp1.nc
    cdo -select,timestep=${ntime1}/${ntime2} $climfile tmp2.nc

    cdo -add tmp2.nc tmp1.nc $fcfile # order important, inherits year from first 
    
    rm tmp1.nc tmp2.nc
  fi

done
