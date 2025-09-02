#!/bin/bash
#
# Usage: ./extract_mom6nep_selected.sh x1 x2 y1 y2 yrstr yrend archdir simname ddmm
#
# This script extracts 

#--------------------
# Setup
#--------------------

# Hyperslab indices (Alaska region)

x1=$1
x2=$2
y1=$3
y2=$4

# Year range to extract

yrstr=$5
yrend=$6

# Where is the MOM6-NEP output?

# arch_dir=/archive/e1n/fre/cefi/NEP/2024_11/NEP_nudge_spinup/gfdl.ncrc6-intel23-repro/history
arch_dir=$7


# Where will the final files be staged?
# Note: Keeping structure similar to old ROMS processing:
#   Level 1 = direct output from model, just rearranged for easier access
#   Level 2 = extra diagnostic variables calculated during postprocessing
#   Level 3 = more highly processed, ususually involving spatially-condensed 
#             indices (e.g. regional average, survey rep)

simname=$8
# simname=mom6nep_hc202411
transferfol=/collab1/data_untrusted/Kelly.Kearney/${simname}/Level1-2/ 

# Variables from ocean, cobalt, and ice files

varo="tos,tob" 
varc="btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf"
vari="siconc"

if [ "$#" -lt 9 ]; then
    ddmm="0101"
else
    ddmm=$9
fi

#--------------------
# Extract
#--------------------

for (( yr=$yrstr; yr<=$yrend; yr++ )); do

    echo $yr

    # Extract specific netcdf files from the tar, extract selected variables and region, delete original

    echo "   extracting ocean variables..."
    # Ocean daily

    tar -xf $arch_dir/${yr}${ddmm}.nc.tar ./${yr}${ddmm}.ocean_daily.nc
    ncks -F -O -d ih,${x1},${x2} -d jh,${y1},${y2} -v ${varo} ${yr}${ddmm}.ocean_daily.nc ${simname}_selected_daily_${yr}.nc
    rm ${yr}${ddmm}.ocean_daily.nc

    # COBALT daily

    echo "   extracting COBALT variables..."
    tar -xf $arch_dir/${yr}${ddmm}.nc.tar ./${yr}${ddmm}.ocean_cobalt_daily_2d.nc
    ncks -F -A -d ih,${x1},${x2} -d jh,${y1},${y2} -v ${varc} ${yr}${ddmm}.ocean_cobalt_daily_2d.nc ${simname}_selected_daily_${yr}.nc
    rm ${yr}${ddmm}.ocean_cobalt_daily_2d.nc

    # Ice daily

    echo "   extracting ice variables..."
    tar -xf $arch_dir/${yr}${ddmm}.nc.tar ./${yr}${ddmm}.ice_daily.nc
    ncrename -d xT,ih -d yT,jh ./${yr}${ddmm}.ice_daily.nc # renaming ice dimensions for easier processing later
    ncks -F -A -d ih,${x1},${x2} -d jh,${y1},${y2} -v ${vari} ${yr}${ddmm}.ice_daily.nc ${simname}_selected_daily_${yr}.nc
    rm ${yr}${ddmm}.ice_daily.nc

    # Move to globus untrusted endopoint staging area
    
    echo "   moving to staging folder..."

    mv ${simname}_selected_daily_${yr}${ddmm}.nc $transferfol

done
