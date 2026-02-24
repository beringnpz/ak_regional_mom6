#!/bin/bash

USEAGE="Usage: ./subset_from_archive.sh [--subdomain <iq1 iq2 jq1 jq2>]
        [--years <yrstr yrend>] [--mmdd <mmdd>] [--archdir <arch_dir>]
        [ftype1 <varstr1> ...] [-h]  

where

  --subdomain <iq1 iq2 jq1 jq2>: subdomain range, with grid variables 
        extracted in the horizontal subregion defined by iq=iq1:iq2, 
        jq = jq11:jq2 (ih = iq1+0.5:iq2-0.5, jh=iq1+0.5:iq2-0.5)
        The specified subdomain must encompass at least one h-point,
        i.e., iq2-iq1>0 and jq2-jq1>0.

  --years <yrstr yrend>: years to process, indicated by starting and ending 
        years of the desired range.  The script expects to find yearly-chunked 
        output files identified by start date (YYYYMMDD); data will be 
        identified by this naming scheme rather than by file contents.
    
  --mmdd <mmdd>: month and day associated with each archive file start date.  
        Default of 0101 applies to most yearly-chunked data.

  --archfdir <arch_dir>: path to archive folder where tarred dataset is found
        (without trailing slash)

  --ppdir <ppfol>: path to processed data folder where output files will be 
        placed.  This will be the current directory unless otherwise specified.

  --ftypeX varstrX: file type (corresponding to a file_name in the MOM6 diag 
        table) and comma-delimited list of variables to extract from each file 
        type.  For example 'ocean_daily tob,tos' will extract tob and tos 
        variables from the MMDDYYYY.ocean_daily.nc files.  Any number of 
        ftype/varstr pairs can be passed as input.

  --region <region>: regional abbreviation that applies to the extracted data, 
        used for file naming only

  --release <release>: release abbreviation (or other preferred simulation name)
        that applies to the extracted data, used for file naming only

This function extracts a subset of variables across the indicated horizontal 
subregion from a MOM6 simulation archive. The resulting files follow the 
following naming scheme, which loosely mimics the CEFI Data Portal conventions:

<ppdir>/<ftype>.<region>.<subdomain>.hcast.<release>.YYYYMMDD.nc
"

#--------------------
# Setup
#--------------------

# Defaults

mmdd="0101"
ppfol="."

# Input parsing

ftype=()
varstr=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --subdomain)
      iq1=$2
      iq2=$3
      jq1=$4
      jq2=$5
      shift # past argument
      shift # past value
      shift
      shift
      shift
      ;;
    --years)
      yrstr=$2
      yrend=$3 
      shift # past argument
      shift # past value
      shift
      ;;    
    --mmdd)
      mmdd=$2
      shift # past argument
      shift # past value
      ;;
    --ppdir)
      ppfol=$2
      shift # past argument
      shift # past value
      ;;
    --archdir)
      arch_dir=$2
      shift
      shift
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
    -h|--help) 
       echo "$USEAGE"
       exit
       ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      ftype+=("$1") # save positional arg
      varstr+=("$2")
      shift # past argument
      shift
      ;;
  esac
done

#--------------------
# Extract
#--------------------

# Subgrid details (note: using 0-based default option for ncks)
#
#   q: iq1-iq2, jq1-jq2
#   h: iq1-(iq2-1), jq1-(jq2-1)
#   u: iq1-iq2, jq1-jq2
#   v: iq1-iq2, jq1-jq2
#
# The smallest "block" (jq2-jq1==1 and iq2-iq1==1) includes 1 h-point, 4 q-points,
# 2 u-points, and 2 v-points.  For simplicity, I'm not allowing subsets that
# don't include at least one h-point.

di=$((iq2-iq1))
dj=$((jq2-jq1))

if !((${#iq1[@]})); then # iq1 not set
    subsetflag=0
    subdomainstr="full"
elif  (( di < 1 )) || (( dj < 1 )) ; then
    echo "Subregion must span at least one h-point"
    exit 1
else
    subdomainstr="iq${iq1}-${iq2}jq${jq1}-${jq2}"
    ih2=$((iq2-1))
    jh2=$((jq2-1))
    subsetflag=1
fi

for (( yr=$yrstr; yr<=$yrend; yr++ )); do

    echo $yr

    for (( i=0; i<${#ftype[@]}; i++ )); do

        echo "   extracting variables from ${ftype[$i]}..."

        # Untar from archive

        tar -xf $arch_dir/${yr}${mmdd}.nc.tar ./${yr}${mmdd}.${ftype[$i]}.nc
    
        # Rename dimensions as needed (note: this does not change the values 
        # of any coordinate variables associated with the dimensions!)

        ncrename -d .xT,ih  -d .yT,jh \
                 -d .xTe,iq -d .yTe,jq \
                 -d .xh,ih  -d .yh,jh \
                 -d .xq,iq  -d .yq,jq ./${yr}${mmdd}.${ftype[$i]}.nc

        # Extract relevant variables using the specified regional subset

        dstr=""
        if [ "${subsetflag}" -eq 1 ]; then

            # ncks annoyingly doesn't support the only-if-exists -d option of 
            # ncrename, so we have to check for dimensions

            if ncdump -h ${yr}${mmdd}.${ftype[$i]}.nc | grep -q "iq = "; then
              dstr="${dstr} -d iq,${iq1},${iq2}"
            fi
            if ncdump -h ${yr}${mmdd}.${ftype[$i]}.nc | grep -q "jq = "; then
              dstr="${dstr} -d jq,${jq1},${jq2}"
            fi
            if ncdump -h ${yr}${mmdd}.${ftype[$i]}.nc | grep -q "ih = "; then
              dstr="${dstr} -d ih,${iq1},${ih2}"
            fi
            if ncdump -h ${yr}${mmdd}.${ftype[$i]}.nc | grep -q "jh = "; then
              dstr="${dstr} -d jh,${jq1},${jh2}"
            fi

        fi

        ncks -O ${dstr} \
                -v ${varstr[$i]} \
                ${yr}${mmdd}.${ftype[$i]}.nc \
                ${ppfol}/${region}.${subdomainstr}.hcast.${release}.${yr}${mmdd}.${ftype[$i]}.nc

        # Delete untarred original (cleanup), move new file to 

        rm ${yr}${mmdd}.${ftype[$i]}.nc

    done
done
