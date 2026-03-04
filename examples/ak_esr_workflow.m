% 2025 ESR Workflow, updated
% 
% This is an updated version of the original 2025 ESR workflow, adapted to
% test the newer version of the ak_regional_mom6 code

%% Data processing

%% ... Extract simulation data (run on PPAN)
%
%  First, run ak_esr_peec_extract.sh to get new and/or updated data.
%  Update these variables to match the release details and location of the
%  extracted data.

% Local "CEFI Portal" options

cpopts = {'portalpath', '~/Documents/Data/CEFI/ak_cefiportal', ...
          'release', 'e202507', ...
          'subdomain', [0 342 446 743]};

staticyyyymmdd = '20240101';

%% ... Create masking variables file

[Scs, Mask] = build_ak_masks(...
    'cpopts', [cpopts 'freq', 'static'], ...
    'expandname', true, ....
    'yyyymmdd', staticyyyymmdd, 'setuponly', false, ...
    'maskonly', true, ...
    'newname', 'ak_masks');

%% ... Calculate climatologies, anomalies, and forecasts