#!/bin/bash

USEAGE="Usage: ./subset_from_archive.sh [--subdomain <iq1 iq2 jq1 jq2>]
        [--years <yrstr yrend>] [--mmdd <mmdd>] [--archdir <arch_dir>]
        [ftype1 <varstr1> ...] [--split] [-h]

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

  --archdir <arch_dir>: path to archive folder where tarred dataset is found
        (without trailing slash)

  --ppdir <ppfol>: path to processed data folder where output files will be
        placed (without trailing slash).  This will be the current directory
        unless otherwise specified.

  --ftypeX varstrX: file type (corresponding to a file_name in the MOM6 diag
        table) and comma-delimited list of variables to extract from each file
        type.  For example 'ocean_daily tob,tos' will extract tob and tos
        variables from the MMDDYYYY.ocean_daily.nc files.  Any number of
        ftype/varstr pairs can be passed as input.

  --region <region>: regional abbreviation that applies to the extracted data,
        used for file naming only

  --release <release>: release abbreviation (or other preferred simulation name)
        that applies to the extracted data, used for file naming only

  --exptype <exptype>: experiment type abbreviation that applies to the
        extracted data, used for file naming only. hcast (hindcast) is the
        default if not included

  --coordfile <coordfile>: path to coordinate variable file. If included,
        geolat*/geolon* variables from this file will be appended to the output
        files (note that this will overwrite any existing geolat*/geolon*)
        variables, so do not include if file already contains those)

  --namescheme <namescheme>: naming scheme.  Options are currently:
        portalv1: mimics CEFI Data Portal circa 02/2025, with nested folders,
                  variable-first file names
        kkflat:   flat file structure (no folders), variable-last file names, 
                  same-length frequency abbreviations

  --split: if included, the output will be split into individual files per
        variable

This function extracts a subset of variables across the indicated horizontal
subregion from a MOM6 simulation archive. The resulting files follow the
following naming scheme, which loosely mimics the CEFI Data Portal conventions:

<ppdir>/<region>.<subdomain>.<exptype>.<freq>.<release>.YYYYMMDD.<ftype>.nc
<ppdir>/<region>.<subdomain>.<exptype>.<freq>.<release>.YYYYMMDD.<variable>.nc
"

#--------------------
# Setup
#--------------------

# Defaults

mmdd="0101"
ppfol="."
splitflag=0
exptype="hcast"
namescheme="portalv1"

# CEFI name tables (short-to-long name lookup)

declare -A tbl_region=(
  [nwa]="nep"
  [nep]="northeast_pacific"
  [arc]="arctic"
  [pci]="pacific_islands"
  [glk]="great_lakes"
)

declare -A tbl_subdomain=(
  [full]="full_domain"
)

declare -A tbl_exptype=(
  [hcast]="hindcast"
  [ss_fcast]="seasonal_forecast"
  [ss_fcast_init]="seasonal_forecast_initialization"
  [ss_refcast]="seasonal_reforecast"
  [dc_fcast_init]="decadal_forecast_initialization"
  [dc_fcast]="decadal_forecast"
  [ltm_proj]="long_term_projection"
)

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
    --coordfile)
      coordfile=$2
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
    --split)
      splitflag=1
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

if ((${#coordfile[@]})) && [[ ! -f "${coordfile}" ]]; then
  echo "Specified coordinate file ${coordfile} does not exist."
  exit 1
fi

if [[ ! -d "${ppfol}" ]]; then
  echo "Output folder ${ppfol} does not exist."
  exit 1
fi

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
    subdomainstr_long=${tbl_subdomain["full"]}
elif  (( di < 1 )) || (( dj < 1 )) ; then
    echo "Subregion must span at least one h-point"
    exit 1
else
    subdomainstr="iq${iq1}-${iq2}jq${jq1}-${jq2}"
    subdomainstr_long=${subdomainstr}
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

        # Determine file frequency by parsing average_DT value

        if ncdump -h ${yr}${mmdd}.${ftype[$i]}.nc | grep -q "average_DT"; then

          freqdays=$( ncdump -v average_DT ${yr}${mmdd}.${ftype[$i]}.nc | grep "average_DT =" )

          freqdays=${freqdays/average_DT =/}
          freqdays=${freqdays/;/}

          if [[ "${freqdays}" == *","*  ]]; then # if multiple (contains comma)
            IFS=',' read -ra freqdays <<< "${freqdays}"
          else
            freqdays=( "${freqdays}" )
          fi

          if (( freqdays[0] == 1 )); then
            freq="day"
            freq_long="daily"
          elif (( freqdays[0] >= 28 )) && (( freqdays[0] <= 31 )); then
            freq="mon"
            freq_long="monthly"
          elif (( freqdays[0] >= 365 )) && (( freqdays[0] <= 366 )); then
            freq="ann"
            freq_long="yearly"
          else
            echo "Could not parse file frequency"
            exit 1
          fi
        else
          freq="sta"
          freq_long="static"
        fi

        # Extract requested variables and region to new file
        
        case ${namescheme} in
          portalv1)
            swapflag=1
            newbase=".${region}.${subdomainstr}.${exptype}.${freq_long}.${release}.${yr}${mmdd}" # for splits
            newfile="${ftype[$i]}${newbase}.nc" # for consolidated
            outfilt="*${newbase}.nc"
            outfol="${ppfol}/${tbl_region[${region}]}/${subdomainstr_long}/${tbl_exptype[${exptype}]}/${freq_long}/raw/${release}"
            ;;
          kkflat)
            newfile="${region}.${subdomainstr}.${exptype}.${freq}.${release}.${yr}${mmdd}.${ftype[$i]}.nc"
            newbase="${region}.${subdomainstr}.${exptype}.${freq}.${release}.${yr}${mmdd}." # for splits
            newfile="${newbase}${ftype[$i]}.nc" # for consolidated
            outfilt="${newbase}*.nc"
            outfol="${ppfol}"
            swapflag=0
            ;;
          *)
            echo "Unrecognized naming scheme ${namescheme}"
            exit 1
        esac

        if [[ ! -d "${outfol}" ]]; then
          mkdir -p ${outfol}
        fi
        
        ncks -O ${dstr} \
                -v ${varstr[$i]} \
                ${yr}${mmdd}.${ftype[$i]}.nc \
                ${newfile}

        # (cleanup) delete untarred original

        rm ${yr}${mmdd}.${ftype[$i]}.nc

        # Determine which grid each variable is on and set coordinate attribute
        # (using create mode, so if coordinate attribute already exists, no
        # modifications are made)

        hflag=0
        cflag=0
        uflag=0
        vflag=0
        
        if [[ "${varstr[$i]}" == *","*  ]]; then # if multiple (contains comma)
          IFS=',' read -ra varnames <<< "${varstr[$i]}"
        else
          varnames=( "${varstr[$i]}" )
        fi
        for vv in "${varnames[@]}"; do

          if ncdump -h ${newfile} | grep " ${vv}(" | grep -q "jh, ih)"; then # h-variable
            ncatted -O -a coordinates,${vv},c,c,"geolat geolon" ${newfile}
            hflag=1
          elif ncdump -h ${newfile} | grep " ${vv}(" | grep -q "jq, iq)"; then # c-variable
            ncatted -O -a coordinates,${vv},c,c,"geolat_c geolon_c" ${newfile}
            cflag=1
          elif ncdump -h ${newfile} | grep " ${vv}(" | grep -q "jh, iq)"; then # u-variable
            ncatted -O -a coordinates,${vv},c,c,"geolat_u geolon_u" ${newfile}
            uflag=1
          elif ncdump -h ${newfile} | grep " ${vv}(" | grep -q "jq, ih)"; then # v-variable
            ncatted -O -a coordinates,${vv},c,c,"geolat_v geolon_v" ${newfile}
            vflag=1
          else
            echo "Variable ${vv} did not match an expected grid type"
            exit 1
          fi

        done

        # Append the coordinate variables as needed

        if ((${#coordfile[@]})); then # if coordfile was passed
          if [ "${hflag}" -eq 1 ]; then
            ncks -A -v geolat,geolon ${coordfile} ${newfile}
          fi
          if [ "${cflag}" -eq 1 ]; then
            ncks -A -v geolat_c,geolon_c ${coordfile} ${newfile}
          fi
          if [ "${uflag}" -eq 1 ]; then
            ncks -A -v geolat_u,geolon_u ${coordfile} ${newfile}
          fi
          if [ "${vflag}" -eq 1 ]; then
            ncks -A -v geolat_v,geolon_v ${coordfile} ${newfile}
          fi
        fi

        # Split by variable (if requested) and move to final output folder

        if [ "${splitflag}" -eq 1 ]; then
          if [ "${swapflag}" -eq 1 ]; then
            cdo splitname,swap ${newfile} ${newbase}
            rm ${newfile}
            mv ${outfilt} ${outfol}
          else
            cdo splitname ${newfile} ${newbase}
            rm ${newfile}
            mv ${outfilt} ${outfol}
          fi
        else
          mv ${newfile} ${outfol}
        fi

    done
done
