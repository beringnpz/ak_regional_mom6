#!/bin/bash

arch_dir=/archive/e1n/fre/cefi/NEP/2024_11/NEP_nudge_spinup/gfdl.ncrc6-intel23-repro/history

# for yr in {1993..2025}; do
for yr in {1993..1993}; do
    echo $yr
    
    # Extract specific netcdf files from the tar file, extract the variables you want, then delete original
    
    # Ocean daily variables
    
    tar -xf $arch_dir/${yr}0101.nc.tar ./${yr}0101.ocean_daily.nc
    ncks -F -O -d xh,1,265 -d yh,442,302 -v tos,tob ${yr}0101.ocean_daily.nc nep202411_selected_daily_$yr.nc
    rm ${yr}0101.ocean_daily.nc
    
    # Append COBALT variables
    
    tar -xvf $arch_dir/${yr}0101.nc.tar ./${yr}0101.ocean_cobalt_daily.nc
    ncks -F -A -d xh,1,265 -d yh,442,302 -v btm_o2,btm_co3_sol_arag,btm_htotal ${yr}0101.ocean_cobalt_daily.nc nep202411_selected_daily_$yr.nc
    rm daily_cobalt_btm_$yr.nc
    
    # Place on ftp server 
    
    # scp nep202411_selected_daily_$yr.nc Kelly.Kearney@ftp.gfdl.noaa.gov:/pub/kelly.kearney/
    # wait
    # rm nep202411_selected_daily_$yr.nc
    
done

