#!/bin/bash

USAGE="Usage: calculate_persis_forecast.sh [--ppdir <ppfol>] 
      [--subdomain <subdomainstr>] [--exptype <exptype>] [--release <release>]
      [--vars <varstr>] 

where

  --ppdir <ppfol>: path to processed data folder where I/O files will be
        placed (without trailing slash).  This will be the current directory
        unless otherwise specified.

  --subdomain <subdomainstr>: subdomain name, being either a short name from 
        the CEFI Data Portal options or an iqX-XjqX-X style subdomain as created 
        by subset_from_archive.sh

  --exptype <exptype>: experiment type abbreviation that applies to the
        data to which calculatons will be applied. hcast (hindcast) is the 
        default if not included

  --release <release>: release abbreviation (or other preferred simulation name)
        that applies to the data to which calculatons will be applied

  --vars <varstr>: comma-delimited list of variables to which calculatons will 
        be applied

  --filedatestr <filedatestr>: YYYYMMDD file-start-time string corresponding to 
        the anomaly file on which to base the forecast  

  --tindex <tindex>: time index of time step within the anomlay file on which to 
        base the forecast.  If not included, this will be the last time step in
        the file.

This function calculates a persistence-based forecast by adding the 
anomaly-from-climatogy value from a single time step to the climatological daily
timeseries.  This expects the climatology and anomaly files to already exist in 
the freq=daily, grid=extra folder corresponding to the indicated simulation, and
places the new output files (labeled with the fcpersis_ prefix) in the same 
location.
"

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
    
    cdo -add tmp2.nc tmp1.nc tmp3.nc # add
   
    cdo setyear,${fcdatestr:0:4} tmp3.nc $fcfile # ensure year reflects anom, not clim
    
    rm tmp1.nc tmp2.nc tmp3.nc
  
  fi

done
