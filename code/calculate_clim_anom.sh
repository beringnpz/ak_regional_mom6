#!/bin/bash
#
# Syntax: ./calculate_clim_anom.sh simname datafol
#

# if [ "$#" -ne 2 ]; then
#     echo "Script expects 2 inputs (simname, datafol), exiting"
#     exit 9
# fi

#-----------------
# Setup
#-----------------

source cefi_name_tables.sh

# Input parsing

while [[ $# -gt 0 ]]; do
  case $1 in
    --subdomain)
      subdomainstr=$2
      shift
      shift
      ;;
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
    --release)
      release=$2
      shift
      shift
      ;;
    --exptype)
      exptype=$2
      shift
      shift
      ;;
    --namescheme)
      namescheme=$2
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
  # (TODO: Need to figure out how to dynamically build that glob in a way cdo accepts)

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


# climfile=${simfol}/Level3/${simname}_daily_clim_1993-2022.nc

#globstr="\'${simfol}/Level1-2/${simname}_selected_daily*.nc\'"

# Example of bash filename expansion: cdo mergetime data_{01..12}.nc output.nc

# if [ ! -f $climfile ] ; then
#   echo "Building 1993-2022 climatology"
#   cdo -ydaymean -selyear,1993/2022 -cat '/gscratch/cicoes/GR011846_reem/CEFI_data/mom6nep_hc202507/Level1-2/mom6nep_hc202507_selected_daily*.nc' $climfile
#   #cdo -ydaymean -selyear,1993/2022 -cat $globstr $climfile
# fi

# # Calculate daily anomaly relative to climatology

# dailyfiles=( ${simfol}/Level1-2/${simname}_selected_daily_*.nc )

# for nepdailyfile in "${dailyfiles[@]}"; do
  
#   anomfile=${nepdailyfile/Level1-2/Level3}
#   anomfile=${anomfile/selected_daily/daily_anomaly}

#   if ! test -f $anomfile; then
#     echo "Calculating anomalies: ${anomfile}"
#     cdo ydaysub $nepdailyfile $climfile $anomfile
#   fi

# done

# Example

# for r_number in {0..24} ; do 
#   r_tag='r'$(printf '%02d' "$r_number")
#   echo "$r_tag"
#   cdo cat -apply,shifttime,$(( $r_number * 5 ))year [ dtr_m_ECEarth_PD_s01"$r_tag"_20??.nc ] dtr_m_ECEarth_PD_s01r00-24_2035-2159.nc
# done
