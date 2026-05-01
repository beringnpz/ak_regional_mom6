function V = ak_variable_info
%AK_VARIABLE_INFO Returns table with AK-specific info for MOM6 variables
%
% V = ak_variable_info
%
% This function returns a shortcut table that holds some of the common info
% associated with variables and used for different types of plots and
% analysis.
%
% Output variables:
%
%   V:  table, row names correspond to MOM6-COBALT variable names, and
%       columns include the following info:
%
%       limanom:    limits for anomaly plots
%
%       limmap:     limits for spatial plots in the AK regions
%
%       cmap:       colormap for this variable, primarily designed to work
%                   with the limmap values as color limits
%
%       limsvyavg:  limits for time series plots of survey-region-averaged
%                   indices
%
%       labelname:  string, text name suitable for labeling axes, subplots,
%                   etc.

% Copyright 2025 Kelly Kearney

V = {...
    "tob"              [-3 3]             [-2 15]         cmocean('thermal')      [-2 8]     "Bottom temperature (\circC)"   
    "cpool2p0"         [-3 3]             [-2 3]          cmocean('-dense', 5)    [0 1]      "2\circC cold pool index"
    "cpool0p0"         [-3 3]             [-2 3]          cmocean('-dense', 5)    [0 1]      "0\circC cold pool index"
    "tos"              [-3 3]             [-2 15]         cmocean('thermal')      [-2 15]    "Surface temperature (\circC)"
    "pH"               [-0.3 0.3]         [7.5 8.1]       cmocean('-curl')        [7.8 8.2]  "pH"
    "btm_o2"           [-180 180]         [0 300]         cmocean('oxy')          [140 350]  "Bottom oxygen (umol/kg)"
    "omega"            [-1 1]             [0.5 1.5]       cmocean('delta')        [0.8 1.5]  "Aragontite saturation state"
    "btm_co3_ion"      [NaN NaN]          [NaN NaN]       parula()                [NaN NaN]  "Bottom [CO3-]"
    "btm_co3_sol_arag" [NaN NaN]          [NaN NaN]       parula()                [NaN NaN]  "Aragonite solubility"
    "btm_htotal"       [NaN NaN]          [NaN NaN]       parula()                [NaN NaN]  "Bottom [H+]"
};

V = cell2table(V(:,2:end), ...
    'variablenames', {'limanom', 'limmap', 'cmap', 'limsvyavg', 'labelname'}, ...
    'rownames', string(V(:,1)));
