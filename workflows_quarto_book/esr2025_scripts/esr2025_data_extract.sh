
# Static file, portalv1
# (produces /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202507/ocean_static.nep.iq0-342jq446-743.hcast.static.e202507.20240101.nc)

subset_from_archive.sh --region nep \
                       --subdomain 0 342 446 743 \
                       --years 2024 2024 \
                       --mmdd 0101 \
                       --archdir /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history \
                       --release e202507 \
                       --namescheme portalv1 \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       ocean_static areacello,deptho,sftof,Coriolis,geolon,geolat,geolon_c,geolat_c,geolon_u,geolat_u,geolon_v,geolat_v,wet,wet_c,wet_u,wet_v,dxt,dyt,dxCu,dyCu,dxCv,dyCv,areacello_cu,areacello_cv,areacello_bu

# Primary hindcast data (1993-2025/06/31)

subset_from_archive.sh --region nep \
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

subset_from_archive.sh --region nep \
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