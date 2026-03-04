
# ./calculate_clim_anom.sh mom6nep_hc202507 /gscratch/cicoes/GR011846_reem/CEFI_data
# ./calculate_persis_forecast.sh 20250101 101 mom6nep_hc202507 /gscratch/cicoes/GR011846_reem/CEFI_data

calculate_clim_anom.sh --region nep \
                       --subdomain iq0-342jq446-743 \
                       --release e202507 \
                       --exptype hcast \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --vars tos,tob,btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf,siconc
