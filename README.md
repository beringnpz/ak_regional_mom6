# Alaska Regional MOM6

This repository holds a collection of code and data designed to develop Alaska-specific analysis of the CEFI MOM6-NEP regional model.  This code was initially written to support the regional model contribution to the 2025 Spring Preview of Ecological and Economic Conditions (PEEC) in Seattle, WA.  We expect it to serve as a rough template for regional workflows associated with the Alaska regional CEFI effort. 

## Repository structure

- code: scripts in a variety of languages (R, bash, Matlab).  These are for the most part written to be called from within that subfolder, either from the command line or from a local R/Matlab environment as appropriate (see workflow below).
- supporting_data: datasets not directly deriving from the MOM6-NEP simulation
- simulation_data: holding area for MOM6-NEP data (original and post-processed).  Currently only holds the pointer file indicating the local access point for this output.

## Simulation output organization

We followed a similar file organization scheme as in our previous ROMS workflow, which loosely follows satellite data processing conventions:

- Level 0: the raw model output. 
- Level 1: same data as Level 0, but with small tweaks to add metadata or reorganize the data (e.g., concatenating or splitting data between files).
- Level 2: Additional variables calculated based on Level 1 data that still falls on the native model grid.
- Level 3: Any level of processing that requires spatial and/or temporal modification, such as regridding, resampling, or calculating of reduced-dimension indices.
- Level 4: Anything derived from a Level 3 variable (e.g., bias corrected versions, fish models, etc.)
Data for this workflow are currently stored on the UW HPC klone machine.  Discussions are ongoing to determine the best open access location for this data to be archived.

## Getting started

To use this code, clone or download a copy to a local machine.  To run scripts that access the simulation ouptut, follow instructions in the `simulation_data/data_folder_readme.md` file to create an appropriate pointer to where data can be accessed on a local machine as needed. (This is a temporary workaround for our current lack of a centralized data management option... hopefully will be improved upon soon!

Note: this code is primarily provided to provide an open and centralized location for codevelopment of tools within the Alaska CEFI Descision Support Team.  It is definitely a work in progress and is not expected to function in a fully robust manner yet.