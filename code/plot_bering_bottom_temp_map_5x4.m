function h = plot_bering_bottom_temp_map_5x4(simname, yrcurrent, opt)
%PLOT_BERING_BOTTOM_TEMP_MAP_5X4 Bering Sea ESR figure
%
% h = plot_bering_bottom_temp_map_5x4(simname, yrcurrent, ...)
%
% This function creates a figure showing the simulated July 1 bottom
% temperature for the last 20 years in the eastern Bering Sea shelf region.
% It is designed to mimic a similar figure produced by the GAP program to
% show bottom temperature measured on the bottom trawl survey
% (https://github.com/afsc-gap-products/coldpool).
%
% Input variables:
%
%   simname:    name of simulation, used to locate output data files.  The
%               path will be constructed as <datafol>/<simname>/Level1-2/.
%
%   yrcurrent:  year for which to produce a plot.  The plot will hold a 5x4
%               grid of axes depicting this year (bottom right) and the
%               previous 19 years' worth of data.
%
% Optional input variables (passed as parameter/value pairs, default in []):
%
%   staticname: base name of static file, assumed to follow the naming
%               scheme <datafol>/<simname>/Level1-2/<staticname>.
%               ["<simname>_ocean_static_ak.nc"]
%
%   datafol:    CEFI data folder path.  Default is the path returned by the
%               cefidatafol.m function

% Copyright 2025 Kelly Kearney

%--------------------
% Parse inputs
%--------------------

arguments
    simname {mustBeTextScalar}
    yrcurrent (1,1) {mustBeInteger}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
end

%--------------------
% Setup
%--------------------

% Read relevant static file variables

staticfile = fullfile(opt.datafol, simname, "Level1-2", opt.staticname);
Grd = ncstruct(staticfile, 'geolat', 'geolon', 'geolat_c', 'geolon_c', ...
    'mask_esr_area', 'mask_survey_strata', 'mask_survey_area', 'wet');

gsz = size(Grd.geolat);

% EBS region

isebs = Grd.mask_esr_area == 3;
latlim = minmax(Grd.geolat(isebs), 'expand');
lonlim = minmax(Grd.geolon(isebs), 'expand');

% Years to plot

yrplt = yrcurrent - (19:-1:0);

% Map decor data: borders, survey region

issebs = Grd.mask_survey_area == 98 & Grd.mask_survey_strata <= 62;
[slon,slat] = mask2poly(Grd.geolon_c, Grd.geolat_c, issebs);

[blat, blon] = deal(cell(1,2));
[blat{1}, blon{1}] = borders('alaska');
[blat{2}, blon{2}] = borders('russia');

warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
bp = polyshape(wrapTo360([blon{:}]), [blat{:}]);
warning(warnstate);

% Projection setup

hfig = figure('visible', 'off');
hb = boxworldmap(latlim, lonlim, 'latgrid', 55:5:60, 'longrid', 180:10:200);
m = getm(hb.ax);
close(hfig);

[xc,yc] = projfwd(m, Grd.geolat_c, Grd.geolon_c);
[bx,by] = projfwd(m, bp.Vertices(:,2), bp.Vertices(:,1));
b = polyshape(bx,by);
[sx, sy] = projfwd(m, slat, slon);

%--------------------
% Plot
%--------------------

% Figure and axis grid creation

h = plotgrid('size', [5 4], 'sp', 0.01, 'mar', 0.02, 'mb', 0.1);
h.ax = h.ax';
setpos(h.fig, '# # 6.5in 9in');
set(h.fig, 'color', 'w');

% Plot bottom temperature

for ii = 1:length(yrplt)

    % Box-ed map axis

    axes(h.ax(ii));
    h.b(ii) = boxworldmap(latlim, lonlim, 'latgrid', 55:5:60, 'longrid', 180:10:200);

    % Read and plot bottom temps for each year

    fglob = fullfile(opt.datafol, simname, "Level1-2", simname+"_selected_daily_" + yrplt(ii) + "*.nc");
    fname = dir(fglob);
    nfile = length(fname);
    fname = fullfile({fname.folder}, {fname.name});
 
    if nfile > 0
        t = ncdateread(fname, 'time');
        [dt,imin] = min(abs(datetime(yrplt(ii),7,1) - t));
        if dt > days(5)
            warning('Possible gap: %s is closest time found', t(imin));
        end
        Tmp = ncstruct(fname, 'tob', struct('time', [imin 1 1]));
        h.p(ii) = pcolorpad(xc, yc, padend(Tmp.tob));
        shading flat;
        uistack(h.p(ii), 'bottom');
    end

    % Add land borders and survey region polygon

    plot(h.ax(ii), bx, by, 'color', rgb('gray'));
    plot(h.ax(ii), sx, sy, 'k');
    set([h.b(ii).lblpar; h.b(ii).lblmer], 'visible', 'off');
end

% Label axes by year

labelaxes(h.ax(1:length(yrplt)), compose('%d',yrplt), 'northwestoutsideabove');

% Set colormap

clim = [-2 3];
cmap = cmocean('-dense', 5);
set(h.ax, 'clim', clim, 'colormap', cmap, 'layer', 'top');
h.cb = colorbar(h.ax(end,1), 'south');
setpos(h.cb, '# 0.02 # #');
xlabel(h.cb, 'Bottom temperature (\circC)');

set(h.ax(length(yrplt)+1:end), 'visible', 'off');

% Add title in lower left

str = 'Bering10K ROMS hindcast, extracted on July 1 of each year';
annotation('textbox', [0.0 0 0.5 0.1], 'string', str, 'vert', 'bottom', ...
    'edgecolor', 'none', 'margin', floor(h.fig.Position(3)*0.02));


end

%--------------------
% Subfunctions
%--------------------

function x = padend(x)
%PADEND Add trailing row and column of NaNs to a 2D array

    x = [x nan(size(x,1),1); nan(1,size(x,2)+1)];

end





