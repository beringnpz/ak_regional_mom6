# Alaska Regional MOM6

This repository holds a collection of code and data designed to develop Alaska-specific analysis of the CEFI MOM6-NEP regional model.  This code was initially written to support the regional model contribution to the 2025 Spring Preview of Ecological and Economic Conditions (PEEC) in Seattle, WA.  We expect it to serve as a rough template for regional workflows associated with the Alaska regional CEFI effort. 

## Repository structure

- code: functions/scripts in a variety of languages (R, bash, Matlab).  These utilities assist in the calculations used in the various workflows when called from the appropriate local R/Matlab/terminal environment.  The plan is to eventually package these as proper toolboxes following the conventions of their respective languages.
- supporting_data: datasets not directly deriving from the MOM6-NEP simulation
- simulation_data: (deprecated) holding area for MOM6-NEP data (original and post-processed).

## Documentation

Workflow documentation, tutorials, etc. can be found in the "Alaska MOM6 regional workflows and skill assessment" document (link coming soon!)

## Note

This code is primarily provided to provide an open and centralized location development of tools within the Alaska CEFI Descision Support Team.  It is definitely a work in progress and is not expected to function in a fully robust manner yet.