# near-real-time extension (2026/01/01 - 2026/08/14)

if [ ! -f $fpart ]; then
subset_from_archive.sh --region nep \
                       --subdomain 0 342 446 743 \
                       --years 2026 2026 \
                       --mmdd 0101 \
                       --archdir /arch0/e1n/fre/cefi/NEP/2026_04/NEP10k_202604_physics_bgc_uv_sponge_august_extension/gfdl.ncrc6-intel25-repro/history/ \
                       --release e202604 \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --namescheme portalv1 \
                       --coordfile /work/Kelly.Kearney/ak_cefiportal/northeast_pacific/iq0-342jq446-743/hindcast/static/raw/e202604/ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc \
                       --split \
                       ocean_daily tos,tob \
                       ocean_cobalt_daily_2d btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf \
                       ice_daily siconc
fi