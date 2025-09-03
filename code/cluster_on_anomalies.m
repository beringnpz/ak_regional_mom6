function [h, Cdata] = cluster_on_anomalies(simname, mmdd, opt)

arguments
    simname {mustBeTextScalar}
    mmdd (1,2) {mustBeInteger}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
    opt.yr {mustBeInteger, mustBeVector} =(1993:year(datetime('today')))
    opt.vars {mustBeText} =["tob", "tos"]
    opt.maskvar ="mask_esr_area"
    opt.maskval =[1 3 4];
    opt.cutoff (1,1) {mustBeNumeric} =0.28;
    opt.Cdata =[]
end

if ischar(opt.vars)
    opt.vars = string(opt.vars);
end

%--------------------
% Setup
%--------------------

staticfile = fullfile(opt.datafol, simname, "Level1-2", opt.staticname);
lev3fol = fullfile(opt.datafol, simname, "Level3");

% Read data for clustering

nyr = length(opt.yr);
nv = length(opt.vars);

%--------------------
% Read data for 
% clustering
%--------------------

% Cluster variable data

if isempty(opt.Cdata)

    for iy = nyr:-1:1
        fprintf('Reading data: %d\n', opt.yr(iy));
    
        % Check for files for this year
    
        fglob = fullfile(opt.datafol, simname, "Level3", simname+"_daily_anomaly_" + opt.yr(iy) + "*.nc");
        fanom = dir(fglob);
        nfile = length(fanom);
        if nfile < 1
            error('No anomaly files found for %d; adjust yr input accordingly', opt.yr(iy));
        end
        fanom = fullfile({fanom.folder}, {fanom.name});
    
        % Find time index corresponding to clustering target date
    
        t = ncdateread(fanom, 'time');
        tidx = interp1(t, 1:length(t), datetime(opt.yr(iy), mmdd(1), mmdd(2)), 'nearest');
    
        Cdata(iy) = ncstruct(fanom, struct('time', [tidx 1 1]), opt.vars{:});
    
    end
else
    Cdata = opt.Cdata;
end

% Mask for regional clustering

Grd = ncstruct(staticfile, 'geolat', 'geolon', 'geolat_c', 'geolon_c', 'wet', opt.maskvar{:});

mask = Grd.wet == 1 & ismember(Grd.(opt.maskvar), opt.maskval);

latlim = minmax(Grd.geolat);
lonlim = minmax(Grd.geolon);

%--------------------
% Cluster
%--------------------

% Reformat cluster variable data into nyr x (nmaskedcell*nprop) array

tmpdata = arrayfun(@(X) cell2mat(cellfun(@(x) x(mask), struct2cell(X), 'uni', 0)), Cdata, 'uni', 0);
tmpdata = cat(2, tmpdata{:})';

% Hierarchical clustering
    
z = linkage(tmpdata, 'ward');

ctmp = clusterdata(tmpdata, ...
    'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
    'distance', 'cutoff', opt.cutoff * max(z(:,3)));

%--------------------
% Plot
%--------------------

% Axis setup

h = plotgrid('size', [nyr 1], 'mr', 0.8, 'sp', 0, 'mar', 0.01);
h.ax = subgridaxes(h.ax, 1, 2);
h.ax = reshape(permute(h.ax, [3 4 2 1]), [], 2);

h.axd = subaxis(1,1,1,'mar', 0.01, 'ml', 0.25);
set(h.fig, 'color', 'w');
setpos(h.fig, '# # 600 1000');

h.axm = axes('position', [h.axd.Position(1)+h.axd.Position(3)-0.25, ...
                          h.axd.Position(2), ...
                          0.25 0.1]);

% Plot map

S = akshapes('dataset', "esr_regions");
Esr = geotable2struct(S.esr_regions(1:4,:));

axes(h.axm);
worldmap(latlim, lonlim);
plotmap(Grd, double(mask));
plotm([Esr.Lat], [Esr.Lon], 'k');
set(h.axm, 'colormap', rgb(["white";"yellow"]), 'clim', [0 1]);
setm(h.axm, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off');

% Plot dendrogram

axes(h.axd)

[h.dh,leafnum,outperm] = dendrogram(z,0, 'label', compose('%d',opt.yr), ...
        'orientation', 'right', 'colorthreshold', opt.cutoff * max(z(:,3)));

set(h.axd, 'ylim', [0.5 nyr+0.5]);

outperm = outperm(end:-1:1);
[~,loc] = ismember(opt.yr, opt.yr(outperm));

% Plot maps

ltlim = minmax(Grd.geolat(mask));
lnlim = minmax(Grd.geolon(mask));

for iy = 1:nyr
    axes(h.ax(loc(iy),1));
    worldmap(ltlim, lnlim);

    tmp = Cdata(iy).tob;
    tmp(~mask) = NaN;
    plotmap(Grd, tmp);

    axes(h.ax(loc(iy),2));
    worldmap(ltlim, lnlim);

    tmp = Cdata(iy).tos;
    tmp(~mask) = NaN;
    plotmap(Grd, tmp);

end

labelaxes(h.ax(loc(end),1), "*", 'northwestoutside');
labelaxes(h.ax(1,:), ["Bottom", "Surface"], 'northoutside', 'fontsize', 6);

arrayfun(@(x) setm(x, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off'), h.ax(1:nyr,:));
set(h.ax, 'clim', [-3 3], 'colormap', cmocean('balance'));
set(h.axd, 'xcolor', 'none');

% labelaxes(h.axd, string(lbl), 'northeast', 'vbuffer', 0);
uistack(h.axm, 'top');

h.cb = colorbar(h.ax(1), 'north');
h.cb.Position = [h.axm.Position(1) h.axm.Position(2)+h.axm.Position(4) h.axm.Position(3) 0.01];
xlabel(h.cb, 'Temperature Anomaly ({\circ}C)');
h.cb.AxisLocation = 'in';

end

function h = plotmap(Grd, var, axflag)
    if nargin < 3
        axflag = false;
    end
    if axflag
        worldmap(minmax(Grd.geolat_c), minmax(Grd.geolon_c));
    end
    h = pcolorm(Grd.geolat_c, Grd.geolon_c, padarray(var, [1 1], NaN, 'post'));
end

function a = geotable2struct(b)
%VERTEXDATA Convert geospatial table to a geographic data structure

    a = table2struct(b);

    for ii = 1:length(a)
        x = a(ii).Shape.InternalData.VertexCoordinate1;
        y = a(ii).Shape.InternalData.VertexCoordinate2;
        idx = a(ii).Shape.InternalData.IndexOfLastVertex;

        n = idx - [1 idx(1:end-1)+1]+1;

        x = mat2cell(x, 1, n);
        y = mat2cell(y, 1, n);
        x = cellfun(@(z) [z NaN], x, 'uni', 0);
        y = cellfun(@(z) [z NaN], y, 'uni', 0);
        x = cat(2, x{:});
        y = cat(2, y{:});

        if isa(a(ii).Shape, 'geopolyshape')
            a(ii).Lat = x;
            a(ii).Lon = wrapTo360(y);
        elseif isa(a(ii).Shape, 'mappolyshape')
            a(ii).Y = y;
            a(ii).X = x;
        end
    end

end

