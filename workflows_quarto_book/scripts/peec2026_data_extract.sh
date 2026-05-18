# TODO: I should rework subset_from_archive.sh so it checks for existence of 
# expected output before running anything.  Until then, I'll do a rough 
# check here.

fstat=/work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202604/ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc
fmain=/work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/daily/raw/e202604/tob.nep.iq0-342jq446-743.hcast.daily.e202604.20250101.nc
fpart=/work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/daily/raw/e202604/tob.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc

# Static file
# Note: should be identical to previous hindcast, but extracting for consistency just in case
# (produces /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202604/ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc

if [ ! -f $fstat ]; then
    subset_from_archive.sh --region nep \
                        --subdomain 0 342 446 743 \
                        --years 2024 2024 \
                        --mmdd 0101 \
                        --archdir /archive/e1n/fre/cefi/NEP/2026_04/NEP10k_202604_physics_bgc_uv_sponge/gfdl.ncrc6-intel25-repro/history/ \
                        --release e202604 \
                        --namescheme portalv1 \
                        --ppdir /work/Kelly.Kearney/ak_cefiportal \
                        ocean_static areacello,deptho,sftof,Coriolis,geolon,geolat,geolon_c,geolat_c,geolon_u,geolat_u,geolon_v,geolat_v,wet,wet_c,wet_u,wet_v,dxt,dyt,dxCu,dyCu,dxCv,dyCv,areacello_cu,areacello_cv,areacello_bu
fi

# Primary hindcast data (1993-2025)

if [ ! -f $fmain ]; then 
    subset_from_archive.sh --region nep \
                        --subdomain 0 342 446 743 \
                        --years 1993 2025 \
                        --mmdd 0101 \
                        --archdir /archive/e1n/fre/cefi/NEP/2026_04/NEP10k_202604_physics_bgc_uv_sponge/gfdl.ncrc6-intel25-repro/history/ \
                        --release e202604 \
                        --ppdir /work/Kelly.Kearney/ak_cefiportal \
                        --namescheme portalv1 \
                        --coordfile /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202604/ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc \
                        --split \
                        ocean_daily tos,tob \
                        ocean_cobalt_daily_2d btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf \
                        ice_daily siconc
fi

# Partial year/Near-real-time (2026) 

if [ ! -f $fpart ]; then
subset_from_archive.sh --region nep \
                       --subdomain 0 342 446 743 \
                       --years 2026 2026 \
                       --mmdd 0101 \
                       --archdir /arch0/e1n/fre/cefi/NEP/2026_04/NEP10k_202604_physics_bgc_uv_sponge_extension/gfdl.ncrc6-intel25-repro/history/ \
                       --release e202604 \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --namescheme portalv1 \
                       --coordfile /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202604/ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc \
                       --split \
                       ocean_daily tos,tob \
                       ocean_cobalt_daily_2d btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf \
                       ice_daily siconc
fi

# Calculate climatology and anomalies

calculate_clim_anom.sh --region nep \
                       --subdomain iq0-342jq446-743 \
                       --release e202604 \
                       --exptype hcast \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --vars tos,tob,btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf,siconc