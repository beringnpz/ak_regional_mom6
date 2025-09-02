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
./extract_mom6nep_selected.sh 1 342 447 743 2025 2025 /archive/e1n/fre/cefi/NEP/2025_07/NEP10k_202507_physics_bgc/gfdl.ncrc6-intel23-repro/history mom6nep_hc202507 0701

# Run on UW Hyak-klone:

# ./calculate_clim_anom.sh mom6nep_hc202507 /gscratch/cicoes/GR011846_reem/CEFI_data/ 2025

