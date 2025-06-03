%% Data analysis and visualization for PEEC

%% Setup

% Data location

datafolfile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'simulation_data', 'data_folder.txt');
if ~exist(datafolfile, 'file')
    error('No data_folder.txt file found');
end
datafol = fileread(datafolfile);

simname = 'mom6nep_hc202411';
lev0fol = fullfile(datafol, simname, 'Level0');
lev1fol = fullfile(datafol, simname, 'Level1-2');
lev3fol = fullfile(datafol, simname, 'Level3');

% The grid (i.e. static) file

grdfileak = fullfile(lev1fol, 'ocean_static_ak.nc');
Grd = ncstruct(grdfileak); % AK only hyperslab
GrdFull = ncstruct(fullfile(lev0fol, '19930101.ocean_static.nc')); % Full region

% Bounds

latlim = minmax(Grd.geolat(~isnan(Grd.mask_esrreg)));
lonlim = minmax(Grd.geolon(~isnan(Grd.mask_esrreg)));

% Borders

bd = load('borderdata.mat');
[~,loc] = ismember({'alaska', 'russia', 'canada'}, lower(bd.places));
blat = [bd.lat{loc}];
blon = [bd.lon{loc}];

% Slide figure setup

outfol = '~/Documents/Conferences/202505_SpringPEEC/slides';
if ~exist(outfol, 'dir')
    mkdir(outfol);
end
savefigs = false; % Flag to export to file

% Silence a few warnings

warning('off', 'map:projections:notStandardProjection');
warning('off', 'MATLAB:polyshape:repairedBySimplify');

%% Slide 1: MOM6 full-grid map
% (Note: slides have 10" x 4.5" usable space, 2.2:1 aspect)
% 
% The first slide depicts the location of the MOM-NEP domain relative to Alaska 
% management regions.

[S, pcrs] = akshapes;

% Map setup

h = plotgrid('size', [1 1], 'mar', 0);
h.fig.Color = 'w';
h.fig.Position(3:4) = [600 600];
plotmap(GrdFull, GrdFull.deptho, true);

proj = getm(gca);

% Transform shapefile polygons to polyshape objects for better plotting

Sps.Esr    = arrayfun(@(x) shape2poly(x, pcrs, proj), S.Esr);
Sps.Survey = arrayfun(@(x) shape2poly(x, pcrs, proj), S.Survey);

% Plot ESR (black) and survey (red) polygons

h.esr = plot(Sps.Esr, 'facecolor', 'w');
h.svy = plot(Sps.Survey(1:end-1), 'facecolor', 'none', 'edgecolor', 'r');

% Prettify

set(gca, 'colormap', cmocean('deep'), 'clim', [0 8000]);
setm(gca, 'frame', 'off', 'plabelmeridian', -120, 'mlabelparall', 22, 'fontsize', 8);
bordersm('countries', 'color', rgb('gray'));

% Add colorbar

h.cb = colorbar('north');
h.cb.Position(1) = 0.7;
h.cb.Position(3) = 0.2;
set(h.cb, 'tickdir', 'out', 'fontsize', 8);
xlabel(h.cb, 'Depth (m)');

if savefigs
    exportgraphics(h.fig, fullfile(outfol, 'mom6nep_map.png'), 'resolution', 150);
end

%% Slide 2: ACC for each variable
%
% Because this is a new simulation, we ran a quick skill assessment to see how 
% useful the persistence-based forecast would be for our various target variables.
% 
%% ... calculate ACC
%
% We calculated anomaly correlation coefficients (ACC) for our target variables, 
% comparing a mid-April persistence-based forecast to the true hindcast value 
% for each year in the 1993-2022 period.
% ACC calculations

% 30-year period for skill assessment

accyr = 1993:2022;

% Target variables

vars = ["tob", "tos", "btm_htotal", "btm_o2"];
vlong = ["Bottom temp. (\circC)", "SST (\circC)", "Bottom pH", "Bottom O_2 (mol/kg)"];

nyr = length(accyr);
nv = length(vars);

% Extract anomalies from Apr 15, June 1, July 1, and Aug 1 of each year

for iy = 1:nyr
    fprintf('%d\n', accyr(iy));
    fanom = fullfile(lev3fol, sprintf('%s_daily_anomaly_%d.nc',simname, accyr(iy)));

    t = ncdateread(fanom, 'time');
    tidx = interp1(t, 1:length(t), datetime(accyr(iy), [4 6 7 8], [15 1 1 1]), 'nearest');
   
    for it = length(tidx):-1:1
        Tmp = ncstruct(fanom, struct('time', [tidx(it) 1 1]));

        for iv = 1:nv
            if iy == 1 && it == length(tidx)
                A.(vars{iv}) = nan([size(Tmp.(vars{iv})) nyr 4]);
            end
            A.(vars{iv})(:,:,iy,it) = Tmp.(vars{iv});
        end
    end
end

% ACC for Apr 15 vs Jun/Jul/Aug 1 (i.e. skill of mid-Apr persistence for
% these 3 summer dates)

for iv = 1:nv
    for it = 1:3
        x = A.(vars{iv});
        Acc.(vars{iv}) = mean(x(:,:,:,1) .* x(:,:,:,2:4),3)./ ...
                         sqrt(mean(x(:,:,:,1).^2,3).*mean(x(:,:,:,2:4).^2,3));
    end
end

%% ... Polyshapes for regional plots
%
% Recalculate polyshapes for the ESR and survey regions using the map projection 
% we'll use for all our AK-only plots.

ftmp = figure('visible', 'off');
atmp = worldmap(latlim, lonlim);
proj = getm(atmp);
Sps.Esr    = arrayfun(@(x) shape2poly(x, pcrs, proj), S.Esr);
Sps.Survey = arrayfun(@(x) shape2poly(x, pcrs, proj), S.Survey);
close(ftmp);

%% ... Plot ACC, one figure per target variable

cmapacc = [cmocean('gray'); crameri('-batlow')];

iplt = 3;

h = plotgrid('size', [2 2], 'sp', -0.01, 'sh', -0.05, 'mar', 0.0, 'ml', 0.3);
h.fig.Position(3:4) = [1300 600];
h.fig.Color = 'w';

for iv = 1:nv
    axes(h.ax(iv));
    worldmap(latlim, lonlim);
    
    vplt = Acc.(vars{iv})(:,:,:,2);
    vplt = {vplt vplt vplt};
    vplt{2}(isnan(Grd.mask_esrreg)) = NaN;
    vplt{3}(isnan(Grd.mask_survey)) = NaN;

    plotmap(Grd, vplt{iplt}, false);
    plotm(blat, blon, 'k');
end

arrayfun(@(x) setm(x, 'grid', 'off', 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off'), h.ax);

set(h.ax, 'clim', [-1 1], 'colormap', cmapacc);

labelaxes(h.ax, vlong, 'northwest', 'fontsize', 16, 'hbuffer', 0);

h.cb = colorbar('south');
h.cb.Position(1) = 0.05;
h.cb.Position(3) = 0.15;
set(h.cb, 'XLim', [-0.25 1], 'fontsize', 14, 'tickdir', 'out');
xlabel(h.cb, 'Anomaly Correlation Coefficient (ACC)', 'fontsize', 16);

if savefigs
    export_fig(h.fig, fullfile(outfol, sprintf('acc_mar15-Jul1_%d',iplt)), '-r150', '-png', '-nocrop');
end

%% Slide 3: Anomaly maps and indices
%
% These slides depict the regionally-averaged timeseries over the entire
% simulation period, including the hindcast period as well as the
% persistence forecast for the current year.  They also display a map of
% the anomaly at the last hindcast time step (mid-April 2025) and the
% persistence forecast that overlays that anomaly on July 1 climatology.

%% ... Prep data for plots

% Read current anomaly and Jul 1 forecast value for all target variables, plus 
% full timeseries for regional indices

% Target variables

vars = ["tob", "tos", "btm_htotal", "btm_o2", "cpool0p0", "cpool2p0", "btm_co3_ion", "btm_co3_sol_arag"];
vlong = ["Bottom temp. (\circC)", "SST (\circC)", "Bottom pH", "Bottom O_2 (mmol/kg)"];

% Forecast and 2025 anomaly data

fcfile = fullfile(lev3fol, 'mom6nep_hc202411_forecast_2025.nc');
anomfile = fullfile(lev3fol, 'mom6nep_hc202411_daily_anomaly_2025.nc');
tfc = ncdateread(fcfile, 'time');
[~,ifc] = min(abs(datetime(2022,7,1)-tfc)); % Note: 2022 year should 2025 in file, update if I fix this

nt = ncinfo(anomfile, 'time');
nt = nt.Size;

Fc = ncstruct(fcfile, struct('time', [ifc 1 1]));
An = ncstruct(anomfile, struct('time', [nt 1 1]));

% Regional index data

idxfile = fullfile(lev3fol, "mom6nep_hc202411_surveyregionavg_" + string(1993:2025)'+".nc");
I = ncinfo(idxfile{1});
idxfcfile = fullfile(lev3fol, 'mom6nep_hc202411_surveyregionavg_forecast_2025.nc');

Tmp = arrayfun(@(x) ncstruct(x, vars{:}, 'time'), idxfile);
for iv = 1:length(vars)
    Idx.(vars{iv}) = cat(1, Tmp.(vars{iv}));
end
Idx.time = cat(1, Tmp.time);
Idx.t = ncdateread(idxfile{1}, 'time', Idx.time);
IdxFc = ncstruct(idxfcfile);
IdxFc.t = ncdateread(idxfcfile, 'time', IdxFc.time);
dv = datevec(IdxFc.t);
dv(:,1) = 2025;
IdxFc.t = datetime(dv);

% Add some additional derived variables (pH, cold pool indices, omega).
% Also, convert oxygen to mmol/kg.

Fc.pH = -log10(Fc.btm_htotal./1.25);
An.pH = Fc.pH + log10((Fc.btm_htotal./1.25 - An.btm_htotal./1.25)); %    log10(abs(An.btm_htotal)) .* sign(An.btm_htotal);
Idx.pH = -log10(Idx.btm_htotal./1.25);
IdxFc.pH = -log10(IdxFc.btm_htotal./1.25);

Fc.cpool0p0 = Fc.tob;
An.cpool0p0 = An.tob;
Fc.cpool2p0 = Fc.tob;
An.cpool2p0 = An.tob;

Fc.btm_o2 = Fc.btm_o2*1e6;
An.btm_o2 = An.btm_o2*1e6;
Idx.btm_o2 = Idx.btm_o2*1e6;
IdxFc.btm_o2 = IdxFc.btm_o2*1e6;

Fc.omega = Fc.btm_co3_ion./Fc.btm_co3_sol_arag;
An.omega = (Fc.btm_co3_ion+An.btm_co3_ion)./(Fc.btm_co3_sol_arag+An.btm_co3_sol_arag) - Fc.omega;
Idx.omega = Idx.btm_co3_ion./Idx.btm_co3_sol_arag;
IdxFc.omega = IdxFc.btm_co3_ion./IdxFc.btm_co3_sol_arag;

%% ... Plot 

regnames = {'Aleutian Islands (AI)', 'Southeast Bering Sea (SEBS)', 'Northern Bering Sea (NBS)', 'Gulf of Alaska (GOA)'};

% Plotting info: variable, color limits for anomaly map,
% color limits for forecast map, colormap, and y-limits for timeseries plot

V = {...
    "tob"            [-3 3]             [-2 15]         cmocean('thermal')      [-2 8]
    "cpool2p0"       [-3 3]             [-2 3]          cmocean('-dense', 5)    [0 1]
    "cpool0p0"       [-3 3]             [-2 3]          cmocean('-dense', 5)    [0 1]
    "tos"            [-3 3]             [-2 15]         cmocean('thermal')      [-2 15]
    "pH"             [-0.3 0.3]         [7.5 8.1]       cmocean('-curl')        [7.8 8.2]
    "btm_o2"         [-180 180]         [0 300]         cmocean('oxy')          [140 350]
    "omega"          [-1 1]             [0.5 1.5]       cmocean('delta')        [0.8 1.5]
};
V = cell2table(V, 'variablenames', {'vplt', 'vlima', 'vlim', 'vmap', 'ylim'});

% One plot per variable 

for iv = 1:height(V)

    % Info from table

    vplt = V.vplt{iv};
    vlima = V.vlima(iv,:);
    vlim = V.vlim(iv,:);
    vmap = V.vmap{iv};

    disp(vplt);

    % Figure setup

    h = plotgrid('size', [4 1], 'mar', 0.05, 'mr', 0.4, 'sp', 0.1);
    h2 = plotgrid('size', [2 1], 'figure', h.fig, 'mar', 0.05, 'ml', 0.63);
    h.fig.Position(3:4) = [1300 600];
    h.fig.Color = 'w';
    
    % Plot anomaly map

    axes(h2.ax(1));
    worldmap(latlim, lonlim);
    plotmap(Grd, An.(vplt));
    set(h2.ax(1), 'clim', vlima, 'colormap', cmocean('balance'));
    h.cb(1) = colorbar('north');
    xlabel(h.cb(1), 'Anomaly: Apr 11, 2025', 'fontsize', 16);
    
    % Plot forecast map

    axes(h2.ax(2));
    worldmap(latlim, lonlim);
    plotmap(Grd, Fc.(vplt));
    set(h2.ax(2), 'clim', vlim, 'colormap', vmap);
    h.cb(2) = colorbar('north');
    xlabel(h.cb(2), 'Prediction: July 1, 2025', 'fontsize', 16);
    
    set(h.cb, 'axisloc', 'out', 'tickdir', 'out', 'fontsize', 14);
    arrayfun(@(x) setm(x, 'grid', 'off', 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off'), h2.ax);
    
    % Plot indices

    for ir = 1:4
        h.idx(ir) = plotindex(h.ax(ir), Idx.t, Idx.(vplt)(:,ir), IdxFc.t, IdxFc.(vplt)(:,ir));
        if ir == 1
            h.leg = legendflex(h.idx(ir).lndoy(end-3:end), cellstr(string(h.idx(ir).yrplt(end-3:end))), ...
                'ref', h.idx(ir).ax(2), 'anchor', {'ne','se'}, 'buffer', [0 0], 'nrow', 2, 'box', 'off', 'fontsize', 12);

            h.leg2 = legendflex([h.idx(ir).lndaily h.idx(ir).lnsummer, h.idx(ir).lnfc], ...
                {'Hindcast (daily)', 'Hindcast (July 1)', 'Forecast'}, ...
                'ref', h.idx(ir).ax(1), 'anchor', {'ne','se'}, 'buffer', [0 0], 'nrow', 1, 'box', 'off', 'fontsize', 12, 'xscale', 0.5);


        end
    end
    labelaxes(arrayfun(@(x) x.ax(1), h.idx), regnames, 'northwestoutsideabove', 'fontsize', 16);
    set([h.idx.ax], 'ylim', V.ylim(iv,:));

    if savefigs
        export_fig(h.fig, fullfile(outfol, sprintf('map_index_%s',vplt)), '-r150', '-png', '-nocrop');
    end
end

%% Additional Analysis/Plots

%% ACC: monthly regional average full grid
% An initial pass at replicating the monthly ACC grids, similar to the Kearney 
% et al., 2021 (DOI: 10.1029/2021jc017545) ROMS Bering10K-based sea ice and bottom 
% temperature predictability study.

Tmp = table2timetable(struct2table(Idx));
Mn = retime(Tmp, 'monthly', 'mean');
Mn = Mn(year(Mn.t) < 2025,:);

regnames = {'Aleutian Islands (AI)', 'Southeast Bering Sea (SEBS)', 'Northern Bering Sea (NBS)', 'Gulf of Alaska (GOA)'};
nreg = length(regnames);

vars = ["tob", "tos", "btm_o2"];
nv = length(vars);

ttmp = reshape(Mn.t, 12, []);

nyr = size(ttmp,2)-2;
nlead = 24;

accmn = nan(12,nlead,nv,nreg);

for iv = 1:nv
    for ir = 1:nreg

        y = Mn.(vars{iv})(:,ir);
        y = reshape(y, 12, []);

        yclim = mean(y, 2);
        yanom = y - yclim;

        for ii = 1:12
            yfc = yanom(ii,1:nyr);
            for il = 1:nlead
                idx = sub2ind(size(yanom), ones(1,nyr)*ii,1:nyr) + il;
                yobs = yanom(idx);
               
               accmn(ii,il,iv,ir) = sum(yfc.*yobs)./(sqrt(sum(yfc.^2))*sqrt(sum(yobs.^2)));
            end
        end
    end
end

% Plot

h = plotgrid('size', [nv nreg], 'sp', 0.01, 'collabel', regnames, 'rowlabel', vars);
for iv = 1:nv
    for ir = 1:nreg
        axes(h.ax(iv,ir));
        pcolorpad(accmn(:,:,iv,ir)');
    end
end
cmapacc = [cmocean('gray', 10); crameri('-batlow', 10)];
set(h.ax, 'colormap', cmapacc, 'clim', [-1 1]);

%% Yearly clustering based on April 11 tos, tob data

clusteryr = 1993:2025;
clustervars = ["tob", "tos"];

nyr = length(clusteryr);
nv = length(clustervars);

for iy = nyr:-1:1
    fprintf('%d\n', clusteryr(iy));
    fanom = fullfile(lev3fol, sprintf('%s_daily_anomaly_%d.nc',simname, clusteryr(iy)));

    % fanom = fullfile(fol, 'pp', sprintf('%s_daily_anomaly_%d.nc','nep202411', accyr(iy)));

    t = ncdateread(fanom, 'time');
    tidx = interp1(t, 1:length(t), datetime(clusteryr(iy), 4, 11), 'nearest');
   
    Cdata(iy) = ncstruct(fanom, struct('time', [tidx 1 1]), clustervars{:});

end

%% ... full grid

cutoff = 0.28;

hasdata = Grd.wet==1;

tmpdata = arrayfun(@(x) [x.tob(hasdata); x.tos(hasdata)], Cdata, 'uni', 0);
tmpdata = cat(2, tmpdata{:})';

z = linkage(tmpdata, 'ward');

ctmp = clusterdata(tmpdata, ...
    'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
    'distance', 'cutoff', cutoff * max(z(:,3)));


hfig = figure('color', 'w');
hfig.Position(3:4) = [600 600];

[hh,t,p] = dendrogram(z,0, 'label', compose('%d',clusteryr), ...
    'orientation', 'right', 'colorthreshold', cutoff * max(z(:,3)));

set(gca, 'xcolor', 'none', 'tickdir', 'out', 'fontsize', 14);
set(hh, 'linewidth', 1.5);

%% ... by region, survey

nreg = max(Grd.mask_survey(:))-1;
regnames = {'Aleutian Islands (AI)', 'Southeast Bering Sea (SEBS)', 'Northern Bering Sea (NBS)', 'Gulf of Alaska (GOA)'};

h = plotgrid('size', [1 nreg], 'mar', 0.05);

for ir = 1:nreg

    hasdata = Grd.wet==1  & Grd.mask_survey == ir;

    tmpdata = arrayfun(@(x) [x.tob(hasdata); x.tos(hasdata)], Cdata, 'uni', 0);
    tmpdata = cat(2, tmpdata{:})';

    z = linkage(tmpdata, 'ward');

    ctmp = clusterdata(tmpdata, ...
    'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
    'distance', 'cutoff', cutoff * max(z(:,3)));
    
    axes(h.ax(ir));
    [hh,t,p] = dendrogram(z,0, 'label', compose('%d',clusteryr), ...
        'orientation', 'right', 'colorthreshold', cutoff * max(z(:,3)));

end
labelaxes(h.ax, regnames, 'northoutside');

%% ... by region, ESR

cutoff = 0.2;

nreg = max(Grd.mask_esrreg(:));
regnames = strtrim(split(ncreadatt(grdfileak, 'mask_esrreg', 'flag_meanings'), ","));

h = plotgrid('size', [1 nreg], 'mar', 0.05);

for ir = 1:nreg

    hasdata = Grd.wet==1  & Grd.mask_esrreg == ir;

    tmpdata = arrayfun(@(x) [x.tob(hasdata); x.tos(hasdata)], Cdata, 'uni', 0);
    tmpdata = cat(2, tmpdata{:})';

    z = linkage(tmpdata, 'ward');

    ctmp = clusterdata(tmpdata, ...
    'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
    'distance', 'cutoff', cutoff * max(z(:,3)));
    
    axes(h.ax(ir));
    [hh,t,p] = dendrogram(z,0, 'label', compose('%d',clusteryr), ...
        'orientation', 'right', 'colorthreshold', cutoff * max(z(:,3)));

end
labelaxes(h.ax, regnames, 'northoutside');

%% ... plot anomalies with dendrograms

nreg = max(Grd.mask_esrreg(:));
regnames = strtrim(split(ncreadatt(grdfileak, 'mask_esrreg', 'flag_meanings'), ","));

regorder = [0 2 1 7 6 5 4 3 -1 -2 -3];

cutoff = 0.28;

for ir = 1:length(regorder)
    if regorder(ir) == 0
        mask = Grd.wet==1;
        lbl = 'Alaska ESR bounding box';
    elseif regorder(ir) == -1
        mask = Grd.wet==1 & ismember(Grd.mask_esrreg, [3 4]);
        lbl = 'Bering Sea (NBS+EBS)';
    elseif regorder(ir) == -2
        mask = Grd.wet==1 & ismember(Grd.mask_esrreg, [1 2 7]);
        lbl = 'Aleutian Islands';
    elseif regorder(ir) == -3
        mask = Grd.wet==1 & ismember(Grd.mask_esrreg, [5 6]);
        lbl = 'Gulf of Alaska';
    else
        mask = Grd.wet==1 & Grd.mask_esrreg == regorder(ir);
        lbl = regnames{regorder(ir)};
    end
    
    tmpdata = arrayfun(@(x) [x.tob(mask); x.tos(mask)], Cdata, 'uni', 0);
    tmpdata = cat(2, tmpdata{:})';
    
    z = linkage(tmpdata, 'ward');
    
    ctmp = clusterdata(tmpdata, ...
        'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
        'distance', 'cutoff', cutoff * max(z(:,3)));
    
    % Axis setup
    
    h = plotgrid('size', [33 1], 'mr', 0.8, 'sp', 0, 'mar', 0.01);
    h.ax = subgridaxes(h.ax, 1, 2);
    h.ax = reshape(permute(h.ax, [3 4 2 1]), [], 2);
    
    h.axd = subaxis(1,1,1,'mar', 0.01, 'ml', 0.25);
    set(h.fig, 'color', 'w');
    setpos(h.fig, '# # 600 1000');

    h.axm = axes('position', [h.axd.Position(1)+h.axd.Position(3)-0.25, ...
                              h.axd.Position(2), ...
                              0.25 0.1]);

    % Plot map

    axes(h.axm);
    worldmap(latlim, lonlim);
    plotmap(Grd, double(mask));
    plotm([S.Esr.Lat], [S.Esr.Lon], 'k');
    set(h.axm, 'colormap', rgb(["white";"yellow"]), 'clim', [0 1]);
    setm(h.axm, 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off', 'grid', 'off');

    % Plot dendrogram
    
    axes(h.axd)
    
    [h.dh,leafnum,outperm] = dendrogram(z,0, 'label', compose('%d',clusteryr), ...
            'orientation', 'right', 'colorthreshold', cutoff * max(z(:,3)));
    
    set(h.axd, 'ylim', [0.5 nyr+0.5]);
    
    outperm = outperm(end:-1:1);
    [~,loc] = ismember(clusteryr, clusteryr(outperm));
    
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

    labelaxes(h.axd, string(lbl), 'northeast', 'vbuffer', 0);
    uistack(h.axm, 'top');

    h.cb = colorbar(h.ax(1), 'north');
    h.cb.Position = [h.axm.Position(1) h.axm.Position(2)+h.axm.Position(4) h.axm.Position(3) 0.01];
    xlabel(h.cb, 'Temperature Anomaly ({\circ}C)');
    h.cb.AxisLocation = 'in';


    if savefigs
        export_fig(h.fig, fullfile(outfol, sprintf('cluster_esrreg%d', regorder(ir))), '-png', '-r150', '-nocrop');
    end

end

%% Subfunctions

function h = plotindex(ax, t, y, tfc, yfc)

    % Split axis

    h.ax = subgridaxes(ax, 1, [0.65 0.05 0.3]);
    h.ax = permute(h.ax, [3 2 1]);
    delete(h.ax(2));
    h.ax = h.ax([1 3]);

    % Full timeseries (with summer values highlighted)

    yrlim = minmax(year([t; tfc]));
    tsummer = datetime(yrlim(1):yrlim(end),7,1);
    ysummer = interp1(t, y, tsummer);

    yfcsummer = interp1(tfc, yfc, tsummer);

    yprc = prctile(ysummer(tsummer<datetime(2023,1,1)), [25 50 75]);
    yavg = mean(ysummer(tsummer<datetime(2023,1,1)));
    ystd = std(ysummer(tsummer<datetime(2023,1,1)));

    h.lndaily = plot(h.ax(1), t, y, 'color', rgb('light gray'));
    hold(h.ax(1), 'on');
    h.lnsummer = plot(h.ax(1), tsummer, ysummer, ...
        'marker', 'o', 'color', rgb('green'), 'markerfacecolor', 'w');

    h.yavg = plot(h.ax(1), minmax(t), ones(2,1).*yavg, 'color', rgb('gray'));
    % h.yprc = plot(h.ax(1), minmax(t), ones(2,1).*yprc, 'color', rgb('gray'), 'linestyle', ':');
    h.ystd = plot(h.ax(1), minmax(t), ones(2,1).*(yavg + [-1 1].*ystd), 'color', rgb('gray'), 'linestyle', '--');


    h.lnfc = plot(h.ax(1), tfc, yfc, 'r');
    h.lnfcsummer = plot(h.ax(1), tsummer, yfcsummer, ...
        'marker', 'o', 'color', 'r', 'markerfacecolor', 'w');

    % Day of year plots

    cmap = rgb(["teal";"purple";"medium blue"]);

    [xg, h.yrplt, tmid] = reshapetimeseries(t, y, 'bin', 'date');

    tk = datetime(1993,1:12,1);

    cidx = mod(h.yrplt(end-3:end-1), 3);
    cidx(cidx==0) = 3;

    h.lndoy = plot(doy(tmid), xg, 'color', rgb('light gray'));
    set(h.lndoy(end), 'color', 'r');
    set(h.lndoy(end-3:end-1), {'color'}, num2cell(cmap(cidx,:),2));

    hold(h.ax(2), 'on');
    h.lndoymean = plot(doy(tmid), mean(xg(:,h.yrplt<2023),2), '--k');

    h.lndoyfc = plot(doy(tfc), yfc, '--r');

    set(h.ax, 'box', 'off', 'tickdir', 'out');
    set(h.ax(2), 'xlim', [0 365], 'xtick', doy(tk), 'xticklabel', datestr(tk,'m'));

    linkaxes(h.ax, 'y');

    set([h.lndaily h.lnsummer h.lndoy' h.lndoyfc h.lnfc h.lnfcsummer], 'linewidth', 1.5)

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

function ps = shape2poly(S, pcrs, proj)

    ps = polyshape(S.X, S.Y);
    [lat, lon] = projinv(pcrs, ps.Vertices(:,1), ps.Vertices(:,2));
    [x,y] = projfwd(proj, lat, lon);
    ps.Vertices = [x y];

end