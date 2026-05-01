# CEFI name tables (short-to-long name lookup)

declare -A tbl_region=(
  [nwa]="northwest_atlantic"
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