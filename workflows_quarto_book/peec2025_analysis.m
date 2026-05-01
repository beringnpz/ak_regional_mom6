% peec2025.m

%%
%| label: fig-peec25-mom6gridmap
%| fig-cap: Slide 1-2: Map of MOM6-NEP grid with Alaska management polygons
%| code-fold: true

% Use a static file from the latest release on the Portal to get model grid
 
C = cefiportalopts('region', 'nep', ...
                   'release', 'r20241015', ...
                   'freq', 'monthly'); 

Grd = readcefigridvars(C, {'geolat_c', 'geolon_c', 'deptho'}, 'expandname', false);
[ltlim(1), ltlim(2)] = bounds(Grd.geolat_c(:));
[lnlim(1), lnlim(2)] = bounds(wrapTo360(Grd.geolon_c(:)));

% Plot model grid-cell depths

h.fig = figure('color', 'w');
h.ax = axes('position', [0 0 1 1]);
h.fig.Position(3:4) = [600 600];

padh = @(x) [x nan(size(x,1),1); nan(1,size(x,2)+1)];
worldmap(ltlim, lnlim);
pcolorm(Grd.geolat_c, Grd.geolon_c, padh(Grd.deptho));

proj = getm(gca);

% Add management shapes
% Mix of projected polyshapes vs geostructures is just for plotting
% purposes) 

S = akshapes;
Esr = geotable2struct(S.esr_regions([1 3 4],:));
for ii = 1:length(Esr)
    [x,y] = projfwd(proj, Esr(ii).Lat, Esr(ii).Lon);
    Esr(ii).pshape = polyshape(x,y);
end

Svy = geotable2struct(S.base_layers_survey_area([1:7 11],:));
for ii = 1:length(Svy)
    [Svy(ii).Lat, Svy(ii).Lon] = projinv(projcrs(3338), Svy(ii).X, Svy(ii).Y);
end

h.esr = plot([Esr.pshape], 'facecolor', 'w');
h.svy = plotm([Svy(1:end-1).Lat], [Svy(1:end-1).Lon],  'r');

% Tweak aesthetics, add borders

set(gca, 'colormap', cmocean('deep'), 'clim', [0 8000]);
setm(gca, 'frame', 'off', 'plabelmeridian', -120, 'mlabelparall', 22, 'fontsize', 8);
bordersm('countries', 'color', rgb('gray'));

% Colorbar

h.cb = colorbar('north');
h.cb.Position(1) = 0.7;
h.cb.Position(3) = 0.2;
set(h.cb, 'tickdir', 'out', 'fontsize', 8);
xlabel(h.cb, 'Depth (m)');


%% Anomaly plots and indices slides
%| label: fig-peec25-indices-and-maps
%| fig-cap: From right to left: variable values over time (grey line = daily hindcast values,
%  red line = persistence forecast values, green circles = July 1 values,
%  red circle = July 1 persistence forecast values); variable values versus
%  day of year, highlighting the most recent 4 years with forecast values in the red dashed line; and maps of the (top)
%  spatial anomaly on the last simulated date and (bottom) predicted value
%  based on this anomaly overlaid on July 1 climatological conditions.
%| code-fold: true

C = cefiportalopts('portalpath', '~/Documents/Data/CEFI/ak_cefiportal/', ...
                   'region', 'nep', ...
                   'release', 'e202507', ...
                   'freq', 'daily', ...
                   'yyyymmdd', '20240101', ...
                   'subdomain', 'iq0-342jq446-743'); 


regnum = [98 143 52 47]; % EBS, NBS, AI, GOA
regname = ["EBS", "NBS", "AI", "GOA"];
reglong = ["Southeast Bering Sea", "Northern Bering Sea", "Aleutian Islands", "Gulf of Alaska"];

% Variable setup: which ones do we plot?

vplt = ["tob", "cpool2p0", "cpool0p0", "tos", "btm_o2"]; %, "pH", "omega"];

% Read surveyavg indices from file

Tbl = readindexdata(C.setopts('subdomain', 'aksvyreg'), ...
    'regnum', regnum, ...
    'vars', vplt);

% Remove leap days from climatology

isleap = month(Tbl.Clim.t) == 2 & day(Tbl.Clim.t) == 29;
Tbl.Clim = Tbl.Clim(~isleap,:);

% Remove post-PEEC-2025 hindcast data (note: this is only necessary since
% we're rerunning the PEEC analysis after an additional summer extension
% was added)

Tbl.Hc = Tbl.Hc(Tbl.Hc.t <= datetime(2025,4,11),:);

% Map setup for anomaly/forecast maps

A = akmapprep('cpopts', C, 'maskvar', 'mask_esr_area', 'maskval', [1 3 4], ...
    'latbuffer', 0.01, 'lonbuffer', 0.01);

nreg = length(regnum);
tilesz = [nreg 5];

tnum = @(r,c) sub2ind(tilesz([2 1]), c, r);

V = ak_variable_info;

% Read map slice data
% TODO: need to handle pH, omega anomalies separately

MapAnom = readmom6mapslice(C, datetime(2025,4,11), ...
    'vars', setdiff(vplt, ["cpool2p0", "cpool0p0"]), 'vartype', 'anomaly');
MapFc = readmom6mapslice(C, datetime(2025,7,1), ...
    'vars', setdiff(vplt, ["cpool2p0", "cpool0p0"]), 'vartype', 'fcpersist');
MapAnom.cpool2p0 = MapAnom.tob;
MapAnom.cpool0p0 = MapAnom.tob;
MapFc.cpool2p0 = MapFc.tob;
MapFc.cpool0p0 = MapFc.tob;

for iv = 1:length(vplt)

    % Figure and tiling setup

    h.fig = figure('color', 'w');
    h.fig.Position(3:4) = [1300 600]; % slide size
    h.t = tiledlayout(tilesz(1), tilesz(2));
    h.t.TileSpacing = 'loose';
    h.t.Padding = 'compact';

    % Plot index timeseries

    for ir = 1:length(regnum)
        h.ax(ir,1) = nexttile(tnum(ir,1), [1 2]);
        h.ts(ir) = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ir, ...
            'style', 'date_vs_val_annual_and_daily');

        h.ax(ir,2) = nexttile(tnum(ir,3));
        h.doy(ir) = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ir, ...
            'style', 'doy_vs_val_lines');
    end

    % Prettify timeseries

    set([h.ts.lnsummer h.ts.lnfcsummer], 'markersize', 4, 'linewidth', 1.5);
    set(h.ax, 'ylim', V{vplt{iv}, 'limsvyavg'}, 'fontsize', 8, 'box', 'off');
    labelaxes(h.ax(:,1), compose("%s (%s)", reglong', regname'), 'northwestoutsideabove', 'fontweight', 'b');
    arrayfun(@(x) uistack(x, 'top'), [h.ts.lnsummer])

    % TODO need to add legends

    % Plot maps

    h.max(1) = nexttile(tnum(1,4), [2 2]);
    h.m(1) = plotmom6akmap(A, MapAnom, ...
        'var', vplt{iv}, 'vartype', 'anomaly', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    h.max(2) = nexttile(tnum(3,4), [2 2]);
    h.m(2) = plotmom6akmap(A, MapFc, ...
        'var', vplt{iv}, 'vartype', 'value', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});

    % Prettify and label maps

    set([h.max], 'fontsize', 6, 'box', 'on');
    
    title(h.max(1), {V{vplt{iv}, 'labelname'}, 'Anomaly: Apr 11, 2025'}, 'fontsize', 12);
    title(h.max(2), 'Prediction: Jul 1, 2025', 'fontsize', 12);

    h.m(1).cb.Position(1) = 0.95;
    h.m(2).cb.Position(1) = 0.95;
    set([h.m.cb], 'axislocation', 'out', 'fontsize', 10, 'tickdir', 'out');
end

%%
%| label: fig-peec25-cluster-allESR
%| fig-cap: Cluster analysis based on surface and bottom temperature on Apr 11 across
%| the three Alaska ESR regions.
%| code-fold: true
%| code-eval: false

[hclus, Cdata] = cluster_on_anomalies(...
    [4 11], ...
    'cpopts', C, ...
    'yr', 1993:2025, ...
    'maskval', [1 3 4], ...
    'ndays', 1);