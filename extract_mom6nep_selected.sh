#!/bin/bash

arch_dir=/archive/e1n/fre/cefi/NEP/2024_11/NEP_nudge_spinup/gfdl.ncrc6-intel23-repro/history

for yr in {1995..2025}; do
    echo $yr

    # Extract specific netcdf files from the tar, extract selected variables and region, delete original

    echo "   extracting ocean variables..."
    # Ocean daily

    tar -xf $arch_dir/${yr}0101.nc.tar ./${yr}0101.ocean_daily.nc
    ncks -F -O -d xh,1,265 -d yh,442,743 -v tos,tob ${yr}0101.ocean_daily.nc nep202411_selected_daily_${yr}.nc
    rm ${yr}0101.ocean_daily.nc

    # COBALT daily

    echo "   extracting COBALT variables..."
    tar -xf $arch_dir/${yr}0101.nc.tar ./${yr}0101.ocean_cobalt_daily_2d.nc
    ncks -F -A -d xh,1,265 -d yh,442,743 -v btm_o2,btm_co3_sol_arag,btm_htotal ${yr}0101.ocean_cobalt_daily_2d.nc nep202411_selected_daily_${yr}.nc
    rm ${yr}0101.ocean_cobalt_daily_2d.nc

    # Ice daily

    echo "   extracting ice variables..."
    tar -xf $arch_dir/${yr}0101.nc.tar ./${yr}0101.ice_daily.nc 
    ncks -F -A -d xT,1,265 -d yT,442,743 -v siconc ${yr}0101.ice_daily.nc nep202411_selected_daily_${yr}.nc
    rm ${yr}0101.ice_daily.nc

    # Move to globus untrusted endopoint staging area

    mv nep202411_selected_daily_${yr}.nc /collab1/data_untrusted/Kelly.Kearney/peec2025/

    # Copy to scp ftp server

#    echo "   copying to ftp server"
#    scp nep202411_selected_daily_${yr}.nc Kelly.Kearney@ftp.gfdl.noaa.gov:/pub/kelly.kearney/
#    wait
#    rm nep202411_selected_daily_${yr}.nc


done
