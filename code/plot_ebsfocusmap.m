function h = plot_ebsfocusmap(simname, yrcurrent, opt)
%PLOT_BERING_BOTTOM_TEMP_MAP_5X4 Bering Sea ESR figure
%
% h = plot_ebsfocusmap(simname, yrcurrent, ...)
%
% This function creates an EBS-focused spatial map of values extracted at a
% specified time across years.  It is the generic-ized version of the 5x4
% ESR spatial map of bottom temperature, extended so one can quickly plot
% any variable, either values or anomalies.
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
%
%   yrfilter:   'all':  plot 1970 to yrcurrent. This range reflects the
%                       older ROMS Bering10K hindcast extent for direct
%                       comparison.  
%               'last20': plot last 20 years from yrcurrent.
%               ['all]
%
%   ncol:       number of columns to use in grid layout [6]
%
%   mmdd:       1 x 2 array holding month and day of month when data will
%               be extracted for each year (e.g., [7 1] will extract July 1
%               data for each year). 
%               [7 1]
%
%   subaxprops: cell array of subaxis parameter/value inputs, used to
%               finesse the positions of subpanel axes
%               [{'sp', 0.01, 'mar', 0.02}]
%               
%   axorder:    'column' or 'row', indicating whether yearly data will be
%               plotted in a column-major or row-major layout.
%               ['column']
%
%   yrlabelloc: string, location of year labels (see labelaxes.m for
%               options)
%               ['southwest']
%
%   cbloc:      string location of colorbar (see colorbar.m for options)
%               ['west']

% Copyright 2025 Kelly Kearney

%--------------------
% Parse inputs
%--------------------

arguments
    simname {mustBeTextScalar}
    yrcurrent (1,1) {mustBeInteger}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
    opt.yrfilter {mustBeTextScalar} ="all"
    opt.ncol (1,1) {mustBeInteger} =6
    opt.var {mustBeTextScalar} ='tob'
    opt.vartype {mustBeTextScalar} ='value'
    opt.mmdd (1,2) {mustBeInteger} =[7 1]
    opt.subaxprops ={'sp', 0.01, 'mar', 0.02}
    opt.axorder {mustBeTextScalar} ="column"
    opt.yrlabelloc {mustBeTextScalar} ="southwest"
    opt.cbloc {mustBeTextScalar} ="west"
end

%--------------------
% Setup
%--------------------

A = akmapprep(simname);

% Years to plot

switch opt.yrfilter
    case "all"
        yrplt = 1970:yrcurrent;
    case "last20"
        yrplt = yrcurrent - (19:-1:0);
end
nyr = length(yrplt);

%--------------------
% Plot
%--------------------

% Figure and axis grid creation

% ncol = 6;
nrow = ceil(nyr/opt.ncol);

h = plotgrid('size', [nrow opt.ncol], opt.subaxprops{:});

switch opt.axorder
    case 'column'
    case 'row'
        h.ax = h.ax';
end
% h.ax = h.ax';
h.fig.Position(3:4) = [8.5 11]*72;
set(h.fig, 'color', 'w');

% Plot bottom temperature

for ii = 1:nyr

    % Box-ed map axis

    axes(h.ax(ii));
    h.b(ii) = boxworldmap(A.latlim, A.lonlim, 'latgrid', 55:5:60, 'longrid', 180:10:200);

    % Read and plot bottom temps for each year

    switch opt.vartype
        case 'anomaly'
            fglob = fullfile(opt.datafol, simname, "Level3", simname+"_daily_anomaly_" + yrplt(ii) + "*.nc");
        case 'value'
            fglob = fullfile(opt.datafol, simname, "Level1-2", simname+"_selected_daily_" + yrplt(ii) + "*.nc");
    end
    fname = dir(fglob);
    nfile = length(fname);
    fname = fullfile({fname.folder}, {fname.name});
 
    if nfile > 0
        t = ncdateread(fname, 'time');
        [dt,imin] = min(abs(datetime(yrplt(ii),opt.mmdd(1),opt.mmdd(2)) - t));
        if dt > days(5)
            warning('Possible gap: %s is closest time found', t(imin));
        end
        Tmp = ncstruct(fname, opt.var, struct('time', [imin 1 1]));
        h.p(ii) = pcolorpad(A.xc, A.yc, padend(Tmp.(opt.var)));
        shading flat;
        uistack(h.p(ii), 'bottom');
    end

    % Add land borders and survey region polygon

    plot(h.ax(ii), A.bx, A.by, 'color', rgb('gray'));
    plot(h.ax(ii), A.sx, A.sy, 'k');
    set([h.b(ii).lblpar; h.b(ii).lblmer], 'visible', 'off');
end

% Label axes by year

labelaxes(h.ax(1:length(yrplt)), compose('%d',yrplt), opt.yrlabelloc);

% Set colormap

V = ak_variable_info;

switch opt.vartype
    case 'anomaly'
        set(h.ax, 'clim', V{opt.var, 'limanom'}, 'colormap',cmocean('balance'), 'layer', 'top');
    case 'value'
        set(h.ax, 'clim', V{opt.var, 'limmap'}, 'colormap', V{opt.var, 'cmap'}{1}, 'layer', 'top');
end
h.cb = colorbar(h.ax(end,end), opt.cbloc);

set(h.ax(length(yrplt)+1:end), 'visible', 'off');

end

%--------------------
% Subfunctions
%--------------------

function x = padend(x)
%PADEND Add trailing row and column of NaNs to a 2D array

    x = [x nan(size(x,1),1); nan(1,size(x,2)+1)];

end





