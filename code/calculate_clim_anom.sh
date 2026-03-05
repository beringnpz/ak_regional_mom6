#!/bin/bash

USAGE="Usage: calculate_clim_anom.sh [--ppdir <ppfol>] 
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

This function calculates a daily 1993-2022 climatology timeseries for the 
specified simulation and variables. It also calculates anomaly-from-climatology 
values for all timesteps of the simulation for those variables.  It produces one 
climatology file per variable; anomaly files will follow the same time-chunking 
as the original raw files.  Files are placed in the freq=daily, grid=extra folder 
corresponding to the indicated simulation.
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

#------------------
# Calculate 
#------------------

# I/O folders

basefol="${ppfol}/${tbl_region[${region}]}/${subdomainstr_long}/${tbl_exptype[${exptype}]}"

rawfol="${basefol}/daily/raw/${release}"
outfol="${basefol}/daily/extra/${release}"

if [[ ! -d "${outfol}" ]]; then
  mkdir -p ${outfol}
fi

# Parse variable names

if [[ "${varstr[$i]}" == *","*  ]]; then # if multiple (contains comma)
  IFS=',' read -ra varnames <<< "${varstr[$i]}"
else
  varnames=( "${varstr[$i]}" )
fi


for vv in "${varnames[@]}"; do
  
  # Calculate daily climatology for 30-year period (1993-2022) 

  climfile="${outfol}/clim_${vv}.${region}.${subdomainstr}.${exptype}.daily.${release}.1993-2022.nc"
  globpattern="${rawfol}/${vv}.${region}.${subdomainstr}.${exptype}.daily.${release}.*.nc"

  if [ ! -f $climfile ] ; then
    echo "Building 1993-2022 climatology for ${vv}"
    cdo -ydaymean -selyear,1993/2022 -cat $globpattern $climfile
  fi

  dailyfiles=( ${globpattern} )

  for fname in "${dailyfiles[@]}"; do
    anomfile=${fname/${rawfol}/${outfol}}
    anomfile=${anomfile/${vv}/anom_${vv}}

    if [ ! -f $anomfile ]; then
      echo "Calculating anomalies: ${anomfile}"
      cdo ydaysub $fname $climfile $anomfile
    fi
  done
done
