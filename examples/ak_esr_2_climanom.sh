
# calculate_clim_anom.sh --region nep \
#                        --subdomain iq0-342jq446-743 \
#                        --release e202507 \
#                        --exptype hcast \
#                        --ppdir /work/Kelly.Kearney/ak_cefiportal \
#                        --vars tos,tob,btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf,siconc

calculate_persis_forecast.sh --region nep \
                       --subdomain iq0-342jq446-743 \
                       --release e202507 \
                       --exptype hcast \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --vars tos,tob,btm_o2,btm_co3_sol_arag,btm_htotal,btm_co3_ion,pco2surf,siconc \
                       --filedatestr 20250101 \
                       --tindex 101


calculate_surveyregion_averages.sh --region nep \
                       --subdomain iq0-342jq446-743 \
                       --release e202507 \
                       --exptype hcast \
                       --ppdir /work/Kelly.Kearney/ak_cefiportal \
                       --maskbase ak_masks \
                       --maskdatestr 20240101

