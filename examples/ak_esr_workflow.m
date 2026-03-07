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

C = cefiportalopts('portalpath', '~/Documents/Data/CEFI/ak_cefiportal', ...
                   'release', 'e202507', ...
                   'subdomain', [0 342 446 743], ...
                   'yyyymmdd', '20240101');

yr = 2025;
esrfol = '~/Documents/Manuscripts/2025_ESR/';

%% ... Create masking variables file

[Scs, Mask] = build_ak_masks(...
    'cpopts', C.setopts('freq','static'), ...
    'expandname', true, ....
    'setuponly', true, ...
    'maskonly', true, ...
    'newname', 'ak_masks');

%% ESR Plots

%% ... Map of bottom temperature, past 20 years

fname = fullfile(esrfol, sprintf('btemp_map_%d-%d', [yr-20+1 yr]));

h = plot_bering_bottom_temp_map_5x4(yr, C);

if ~exist([fname '.png'], 'file')
    export_fig(fname, h.fig, '-png', '-r150', '-nocrop');
end