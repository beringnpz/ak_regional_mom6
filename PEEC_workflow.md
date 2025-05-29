# Spring PEEC Workflow

## The simulation

For this simulation, we used a custom extension of the MOM6-NEP hindcast (r202411) simulation.  This simulation was an "unofficial" variant on the primary MOM6-NEP hindcast (accessible on the CEFI portal), and included a few caveats:

- At the time of the extension, the GLORYS dataset that is used for ocean boundary forcing only extended through 03/25/2025, while the atmospheric forcing from ERA5 extended to 04/12/2025, and the GloFAS river forcing reached 04/16/2025.  For the gap between 03/25 and 4/12, the following method was applied for ocean boundary conditions: 
    1) daily climatology for 1995-2024 as constructed
    2) mean delta added from the last 30 days of available GLORYS obcs relative to the same days in the daily climatology, and 
    3) weighting applied to the the last GLORYS day over 10 days into the into the climatology to smooth the transition.
- The simulation originally used a version of the GloFAS-derived river forcing that unintentionally removed the Columbia River.  This bug was fixed for our extension (since the buggy input files were not retained), hence the Columbia "turns on" in 2024.  We do not expect this to affect the Alaska region.
- The boundary sponges used an older algorithm; while since improved upon, the older method shouldn't be detrimental.

This simulation was run by Liz Drenkard at GFDL on 4/22/2025 and archived on GFDL's PPAN (processing and analysis) cluster.

## Analysis

### Initial setup

From my local machine, `export_akgfmaps_polygons.R` was run to export regional shapefiles to the `supporting_data` folder, followed by `build_ocean_static_ak.m` to extract a portion of the MOM6-NEP static file limited to the bounds of Alaska management regions.  To this file Alaska-oriented copy of the static file, we also added a few custom variables to mask the MOM6-NEP h-grid points with AFSC management polygons, including the Ecosystem Status Report (ESR) polygons for the Bering Sea, Aleutian Islands, and Gulf of Alaska, and the groundfish survey regions for these three regions.

### Output extraction (Level 0-1 variables)

`extract_mom6nep_selected.sh` was run on the PPAN system to extract the specific output variables required for this analysis and stage it for Globus transfer to a local analysis machine (processing on PPAN by our team was limited due to a 10GB quota).  The script subsetted the output to the same horizontal hyperslab as for the static file.  It also moved the ice variables from the ice grid to the ocean grid (the grids are identical except in dimension names) to facilitate easier post-processing.

### Calculate derived (Level 3) variables

The `calculate_clim_anom.sh` script calculated the following additional variables, all saved to the Level 3 output folder:

- daily 1993-2022 climatology for target variables
- daily anomaly from climatology values for all target variables 
- a persisence forecast that added the anomaly from the last model time step to the daily climatology for all target varibiables
survey region averages that spatially averaged target variable values across the groundfish survey regions.

### PEEC figures and analysis

The `peec2025_mom6.m` script ran the initial analysis and created figures used in the PEEC slide presentation.  