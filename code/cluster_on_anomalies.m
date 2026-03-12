function [h, Cdata] = cluster_on_anomalies(mmdd, opt)
%CLUSTER_ON_ANALOMALIES Yearly cluster analysis based on spatial anomalies
%
% [h, Cdata] = cluster_on_anomalies(simname, mmdd, opt)
%
% This function performs a simple cluster analysis that identifies
% similarities between years based on the spatial anomalies in each year,
% and creates a dendrogram plot of the results.
%
% This cluster analysis is primarily used to find analogue years from the
% historical record relative to a current management year. While designed
% to support Alaska-region products for the Spring PEEC and summer/fall
% Ecosystem Status Report, this function can be applied more generically,
% assuming MOM6-based simulation output is similarly organized.
%
% Input variables:
%
%   mmdd:       1 x 2 array holding month and day of month when data will
%               be extracted for each year (e.g., [7 1] will extract July 1
%               data for each year).
%
% Optional input variables (passed as parameter/value pairs), default in []:
%
%   cpopts:     cefi portal options object corresponding to the simulation
%
%   ndays:      number of days following the mmdd date over which anomaly
%               data will be averaged for each year.  For example, if mmdd
%               = [7 1] and ndays = 1, the July 1 anomaly will be used.  If
%               ndays=31, the average across July 1-July 31 will be used
%               instead. [1]
%
%   yr:         array of years to include in cluster analysis
%               [1993 through current year]
%
%   vars:       array of strings, anomaly variables to include in cluster
%               analysis ["tob", "tos"]
%
%   maskvar:    masking variable indicating spatial region for cluster
%               analysis.  Can be either a scalar string corresponding to
%               the name of a mask variable in the static file (which will
%               be paired with the maskval input to create a logical array)
%               or a logical array matching the spatial dimensions of the
%               anomaly variables. 
%               ["mask_esr_area"]
%
%   maskval:    Values in the mask array that are included in the analysis.
%               [1 3 4].
%
%   cutoff:     Fractional cutoff value to use for color-coding clusters in
%               the dendrogram display of the cluster analysis.  Clustering
%               will be applied at this fraction times the maximum Ward
%               linkage distance calculated. [0.28]
%
%   Cdata:      cluster data array (see output).  This option is primarily
%               included to allow repeated clustering using the same
%               simulation data but different masking criteria, without
%               needing to reread the data from file.  We note that this
%               option does not come with checks to ensure that this input
%               data matches what would be returned by the specified input
%               variables; please use with caution!    
%
% Output variables:
%
%   h:          structure of handle graphics objects corresponding to the
%               resulting dendrogram plot
%   
%   Cdata:      structure array, with field corresponding to the specified
%               anomaly variable names, each holding a nx x ny x nyr array,
%               where nx and ny are the spatial dimensions of the MOM6
%               regional domain and nyr is the number of years in the
%               analysis.

% Copyright 2025 Kelly Kearney

arguments
    mmdd (1,2) {mustBeInteger}
    opt.cpopts (1,1) {mustBeA(opt.cpopts, "cefiportalopts")} =cefiportalopts()
    opt.ndays =1
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

% Read data for clustering

nyr = length(opt.yr);
nv = length(opt.vars);

%--------------------
% Read data for 
% clustering
%--------------------

% Cluster variable data

if isempty(opt.Cdata)

    for iv = 1:nv
        for iy = nyr:-1:1
            fprintf('Reading data: %s, %d\n', opt.vars{iv}, opt.yr(iy));
    
            % Check for files for this year
                    
            fanom = opt.cpopts.setopts('freq','daily','grid','extra').cefifilelist("anom_"+opt.vars{iv}, opt.yr(iy)+"*");

            % Find time index corresponding to clustering target date
        
            t = ncdateread(fanom, 'time');
    
            tidx = interp1(t, 1:length(t), datetime(opt.yr(iy), mmdd(1), mmdd(2)), 'nearest');
        
            Tmp = ncstruct(fanom, struct('time', [tidx opt.ndays 1]), opt.vars{iv});
            if opt.ndays>1
                Cdata(iy).(opt.vars{iv}) = mean(Tmp.(opt.vars{iv}), 3, 'omitnan');
            else
                Cdata(iy).(opt.vars{iv}) = Tmp.(opt.vars{iv});
            end
        end
    end
else
    Cdata = opt.Cdata;
end

% Mask for regional clustering

if islogical(opt.maskvar) && isequal(size(Cdata(1).(opt.vars{1})), size(opt.maskvar))
    mask = opt.maskvar;
    Grd = readcefigridvars(opt.cpopts, {'geolat', 'geolon', 'geolat_c', 'geolon_c'});
elseif ischar(opt.maskvar) || isstring(opt.maskvar)
    Grd = readcefigridvars(opt.cpopts, {'geolat', 'geolon', 'geolat_c', 'geolon_c', 'wet', opt.maskvar{:}});
    mask = Grd.wet == 1 & ismember(Grd.(opt.maskvar), opt.maskval);
else
    error('Incorrect format for maskvar input')
end


% latlim = minmax(Grd.geolat);
% lonlim = minmax(Grd.geolon);

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

% Anomaly map limits

masklatlim = minmax(Grd.geolat(mask));
masklonlim = minmax(Grd.geolon(mask));

% Figure out anomaly map aspect ratio

hftmp = figure('visible', 'on');
hatmp = axesm('MapProjection', 'eqdconicstd', ...
              'MapLatLimit', masklatlim, ...
              'MapLonLimit', masklonlim);
tightmap;
pos = plotboxpos(hatmp);
aspect = (pos(3)*hftmp.Position(3))./(pos(4).*hftmp.Position(4)); % width/height
close(hftmp);

% Arrange axes: staggered grid of yearly anomaly maps on the left (h.ax),
% dendrogram right (h.axd), inset reference map lower right (h.axm)

figsz = [800 600];

mar = 20; % primary margin, pixels
wdendro = 400; % width dendrogram, pixels
sp1 = 50; % space between maps and dendro, pixels
sp2 = 10; % horizontal space between staggered axes, pixels
scl = 1.5; % scale to increase axes by (>1, <2)
wmasks = figsz(1)-mar*2-sp1-wdendro;

h = plotgrid('size', [nyr 1], ...
             'sv', 0.0, ...
             'mar', mar./figsz(2), ...
             'mt', mar./figsz(2)*2, ...
             'mr', (wdendro+sp1+mar)./figsz(1), ...
             'ml', mar./figsz(1));
h.fig.Position(3:4) = figsz;
set(h.fig, 'color', 'w');

hgt = h.ax(1).Position(4).*figsz(2).*scl; % height, pixels
wdt = hgt.*aspect; % width of each variable sub-panel map, pixels
xcol1 = ((mar+wmasks)-(wdt*nv))./figsz(1);
xcol2 = ((mar+wmasks)-(wdt*nv)*2-sp2)./figsz(1);

for ii = 1:nyr
    if ~mod(ii,2)
        h.ax(ii).Position(1) = xcol2;
    else
        h.ax(ii).Position(1) = xcol1;
    end
    h.ax(ii).Position(3) = wdt*nv./figsz(1);
    h.ax(ii).Position(2) = h.ax(ii).Position(2) + h.ax(ii).Position(4)/2 - hgt/2/figsz(2);
    h.ax(ii).Position(4) = hgt/figsz(2);
end

if nv > 1
    h.ax = subgridaxes(h.ax, 1, nv);
    h.ax = reshape(permute(h.ax, [3 4 2 1]), [], 2);
end

h.axd = subaxis(1,1,1,'mar', mar./figsz(2), 'ml', (wmasks+sp1+mar)./figsz(1), 'mr', mar./figsz(1), 'mt', mar./figsz(2)*2);

h.axm = axes('position', [h.axd.Position(1)+h.axd.Position(3)-0.25, ...
                          h.axd.Position(2), ...
                          0.25 0.1]);

% Add connectors

for ii = 1:nyr
    annotation('line', [h.ax(ii,nv).Position(1)+h.ax(ii,nv).Position(3) (mar+wmasks+sp1/3)./figsz(1)],...
                       (h.ax(ii,nv).Position(2)+h.ax(ii,nv).Position(4)/2).*[1 1], ...
                       'linestyle', ':', 'color', [1 1 1]*0.5);
end

% Plot reference map (locates cluster mask within AK ESR regions)

S = akshapes('dataset', "esr_regions");
Esr = geotable2struct(S.esr_regions(1:4,:));

axes(h.axm);
axesm('MapProjection', 'eqdconicstd', ...
      'MapLatLimit', boundsvec([Esr.Lat]), ...
      'MapLonLimit', boundsvec(wrapTo360([Esr.Lon])));
tightmap;
% worldmap(latlim, lonlim);
plotmap(Grd, double(mask));
plotm([Esr.Lat], [Esr.Lon], 'k');
set(h.axm, 'colormap', rgb(["white";"yellow"]), 'clim', [0 1], 'visible', 'off');
% setm(h.axm, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off');

% Plot dendrogram

axes(h.axd)

[h.dh,leafnum,outperm] = dendrogram(z,0, 'label', compose('%d',opt.yr), ...
        'orientation', 'right', 'colorthreshold', opt.cutoff * max(z(:,3)));

set(h.axd, 'ylim', [0.5 nyr+0.5]);

outperm = outperm(end:-1:1);
[~,loc] = ismember(opt.yr, opt.yr(outperm));

% Change dendrogram colormap from gaudy HSV to more subtle
% (https://tsitsul.in/blog/coloropt/)

xgfs_normal6 = [64, 83, 211; 221, 179, 16; 181, 29, 20; 0, 190, 255; 251, 73, 176; 0, 178, 93; 202, 202, 202];
xgfs_normal12 = [235, 172, 35; 184, 0, 88; 0, 140, 249; 0, 110, 0; 0, 187, 173; 209, 99, 230; 178, 69, 2; 255, 146, 135; 89, 84, 214; 0, 198, 248; 135, 133, 0; 0, 167, 108; 189, 189, 189];
% xgfs_bright6 = [239, 230, 69; 233, 53, 161; 0, 227, 255; 225, 86, 44; 83, 126, 255; 0, 203, 133; 238, 238, 238];
% xgfs_dark6 = [0, 89, 0; 0, 0, 120; 73, 13, 0; 138, 3, 79; 0, 90, 138; 68, 53, 0; 88, 88, 88];
% xgfs_fancy6 = [86, 100, 26; 192, 175, 251; 230, 161, 118; 0, 103, 138; 152, 68, 100; 94, 204, 171; 205, 205, 205];
% xgfs_tarnish6 = [39, 77, 82; 199, 162, 166; 129, 139, 112; 96, 78, 60; 140, 159, 183; 121, 104, 128; 192, 192, 192];

dcol = unique(cat(1, h.dh.Color), 'rows');
dcol = dcol(~ismember(dcol, [0 0 0], 'rows'),:);
ncol = size(dcol,1);

if ncol <= 6
    dcol2 = xgfs_normal6(1:ncol,:);
elseif ncol <= 12
    dcol2 = xgfs_normal12(1:ncol,:);
else
    dcol2 = repmat(xgfs_normal12, ceil(ncol/12),1);
    dcol2 = dcol2(1:ncol,:);
end
dcol2 = dcol2./255;
for ii = 1:ncol
    set(findall(h.dh, 'color', dcol(ii,:)), 'color', dcol2(ii,:));
end

% Plot year-by-year anomaly maps

for iy = 1:nyr
    for iv = 1:nv
        axes(h.ax(loc(iy),iv));
        % worldmap(masklatlim, masklonlim);
        axesm('MapProjection', 'eqdconicstd', ...
              'MapLatLimit', masklatlim, ...
              'MapLonLimit', masklonlim);
        tightmap;

        tmp = Cdata(iy).(opt.vars{iv});
        tmp(~mask) = NaN;
        plotmap(Grd, tmp);
    end
end

V = ak_variable_info;
lbl = strrep(V.labelname(opt.vars), "temperature", "temp.");

lblprop = {'fontsize', 6};
htmp = uicontrol('style', 'text', lblprop{:});
htmp.Position(3:4) = [wdt mar];
lbl = arrayfun(@(x) textwrap(htmp, x), lbl, 'uni', 0);
delete(htmp);

labelaxes(h.ax(1,:), lbl, 'northoutside', lblprop{:});

arrayfun(@(x) setm(x, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off'), h.ax(1:nyr,:));
set(h.ax, 'clim', [-3 3], 'colormap', cmocean('balance'), 'xcolor', [1 1 1]*0.8, 'ycolor', [1 1 1]*0.8);
set(h.ax(loc(end),:), 'xcolor', 'k', 'ycolor', 'k'); % highlight current year
set(h.axd, 'xcolor', 'none');


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

function a = boundsvec(b)
    [a(1), a(2)] = bounds(b);
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

