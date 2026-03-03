# This script provides a record of my MOM6-NEP data extraction for AK regional use
#
# It is intended to be updated regularly as the workflow is updated and extended to
# maintain a rough record of how each tool was called (and to simplify rerunning).

#--------------------
# Spring PEEC 2025
#--------------------

# Run on GFDL PPAN analysis cluster:

# ./extract_mom6nep_selected.sh 1 254 447 743 1993 2024 /archive/e1n/fre/cefi/NEP/2024_11/NEP_nudge_spinup/gfdl.ncrc6-intel23-repro/history mom6nep_hc202411
# Run on UW Hyak-klone:
# (Note: script was originally hard-coded for the mom6nep_hc202411 simulation and klone file structure, and included forecast calculation. 
#  This call is no longer forward-compatible with the newer version)

# ./calculate_clim_anom.sh 2024

#--------------------
# Summer ESR 2025
#--------------------

# Run on GFDL PPAN analysis cluster:

# ./extract_mom6nep_selected.sh 1 342 447 743 1993 2025 /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history mom6nep_hc202507
# ./extract_mom6nep_selected.sh 1 342 447 743 2025 2025 /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history mom6nep_hc202507 0701

# Run on UW Hyak-klone:

# ./calculate_clim_anom.sh mom6nep_hc202507 /gscratch/cicoes/GR011846_reem/CEFI_data
# ./calculate_persis_forecast.sh 20250101 101 mom6nep_hc202507 /gscratch/cicoes/GR011846_reem/CEFI_data

#--------------------
# Updated
# Summer ESR 2025
#--------------------

# Run on GFDL PPAN analysis cluster, using portalv1 naming scheme

# Static file, portalv1
# (produces /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202507/ocean_static.nep.iq0-342jq446-743.hcast.static.e202507.20240101.nc)

./subset_from_archive.sh --region nep \
                         --subdomain 0 342 446 743 \
                         --years 2024 2024 \
                         --mmdd 0101 \
                         --archdir /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history \
                         --release e202507 \
                         --namescheme portalv1 \
                         --ppdir /work/Kelly.Kearney/ak_cefiportal \
                         ocean_static areacello,deptho,sftof,Coriolis,geolon,geolat,geolon_c,geolat_c,geolon_u,geolat_u,geolon_v,geolat_v,wet,wet_c,wet_u,wet_v,dxt,dyt,dxCu,dyCu,dxCv,dyCv,areacello_cu,areacello_cv,areacello_bu

# Primary hindcast data (1993-2025/06/31)

./subset_from_archive.sh --region nep \
                         --subdomain 0 342 446 743 \
                         --years 1993 2025 \
                         --mmdd 0101 \
                         --archdir /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history \
                         --release e202507 \
                         --ppdir /work/Kelly.Kearney/ak_cefiportal \
                         --namescheme portalv1 \
                         --coordfile /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202507/ocean_static.nep.iq0-342jq446-743.hcast.static.e202507.20240101.nc \
                         --split \
                         ocean_daily tos,tob \
                         ocean_cobalt_daily_2d btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf \
                         ice_daily siconc

# near-real-time extension (2025/07/01+)

./subset_from_archive.sh --region nep \
                         --subdomain 0 342 446 743 \
                         --years 2025 2025 \
                         --mmdd 0701 \
                         --archdir /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history \
                         --release e202507 \
                         --ppdir /work/Kelly.Kearney/ak_cefiportal \
                         --namescheme portalv1 \
                         --coordfile /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202507/ocean_static.nep.iq0-342jq446-743.hcast.static.e202507.20240101.nc \
                         --split \
                         ocean_daily tos,tob \
                         ocean_cobalt_daily_2d btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf \
                         ice_daily siconc
