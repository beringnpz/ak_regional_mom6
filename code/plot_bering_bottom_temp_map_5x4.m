function h = plot_bering_bottom_temp_map_5x4(yrcurrent, cpopts)
%PLOT_BERING_BOTTOM_TEMP_MAP_5X4 Bering Sea ESR figure
%
% h = plot_bering_bottom_temp_map_5x4(...)
%
% This function creates a figure showing the simulated July 1 bottom
% temperature for the last 20 years in the eastern Bering Sea shelf region.
% It is designed to mimic a similar figure produced by the GAP program to
% show bottom temperature measured on the bottom trawl survey
% (https://github.com/afsc-gap-products/coldpool).
%
% Input variables:
%
%   yrcurrent:  year for which to produce a plot.  The plot will hold a 5x4
%               grid of axes depicting this year (bottom right) and the
%               previous 19 years' worth of data.
%
% Optional input variables (passed as parameter/value pairs, default in []):
%
%   cpopts:     cefiportalopts object corresponding to the original
%               simulation to be plotted
%               [cefiportalopts()]

% Copyright 2025 Kelly Kearney

%--------------------
% Parse inputs
%--------------------

arguments
    yrcurrent (1,1) {mustBeInteger}
    cpopts (1,1) {mustBeA(cpopts, "cefiportalopts")} =cefiportalopts()
end

%--------------------
% Setup
%--------------------

h = plot_ebsfocusmap('yrcurrent', yrcurrent, ...
    'cpopts', cpopts, ...
    'yrfilter', 'last20', ...
    'ncol', 4, ...
    'var', 'tob', ...
    'vartype', 'value', ...
    'mmdd', [7 1], ...
    'axorder', 'row', ...
    'subaxprops', {'sp', 0.01, 'mar', 0.02, 'mb', 0.1}, ...
    'yrlabelloc', 'northwestoutsideabove', ...
    'cbloc', 'south');

% Aesthetic tweaks for ESR...

h.fig.Position(3:4) = [6.5 9]*72;

% Add title in lower left

str = 'MOM6-COBALT-NEP10k hindcast, extracted on July 1 of each year';
annotation('textbox', [0.0 0 0.5 0.1], 'string', str, 'vert', 'bottom', ...
    'edgecolor', 'none', 'margin', floor(h.fig.Position(3)*0.02));

% Colormap/bar adjustments

clim = [-2 3];
cmap = cmocean('-dense', 5);
set(h.ax, 'clim', clim, 'colormap', cmap);
setpos(h.cb, '# 0.02 # #');
xlabel(h.cb, 'Bottom temperature (\circC)');





