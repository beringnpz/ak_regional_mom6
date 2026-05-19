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
                   'release', 'e202604', ...
                   'freq', 'daily', ...
                   'yyyymmdd', '20240101', ...
                   'subdomain', 'iq0-342jq446-743'); 


regnum = [98 143 52 47]; % EBS, NBS, AI, GOA
regname = ["EBS", "NBS", "AI", "GOA"];
reglong = ["Southeast Bering Sea", "Northern Bering Sea", "Aleutian Islands", "Gulf of Alaska"];

% Variable setup: which ones do we plot?

vplt = ["tob", "cpool2p0", "cpool0p0", "tos", "siconc"]; %, "pH", "omega"];

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

    h.leg(1) = legend([h.ts(1).lndaily, h.ts(1).lnsummer, h.ts(1).lnfc], ...
        ["Hindcast (daily)", "Hindcast (Jul 1)", "Forecast"], ...
        'numcolumns', 2, 'location', 'northeast', 'box', 'off');
    h.leg(1).Position(2) = h.ax(1,1).Position(2)+h.ax(1,1).Position(4);
    h.leg(1).Position(1) = h.ax(1,1).Position(1)+h.ax(1,1).Position(3)-h.leg(1).Position(3);
    
    h.leg(2) = legend(h.doy(1).lndoy(end-3:end), string(2022:2025), ...
        'numcolumns', 2, 'location', 'northwest', 'box', 'off');
    h.leg(2).Position(2) = h.ax(1,2).Position(2)+h.ax(1,2).Position(4);
    h.leg(2).Position(1) = h.ax(1,2).Position(1);


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



%% Animated version

iv = 1; % tob

%-------------------
% Repeat static map
%-------------------

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

%-------------------
% Animate
%-------------------

tplt = datetime(2024,7,1):datetime(2025,7,1);
tplt.Format = 'MMM dd, uuuu';

for ir = 1:length(regnum)
    h.timebar(ir,1) = line(tplt([1 1]), h.ax(ir,1).YLim, 'parent', h.ax(ir,1), 'color', rgb('gold'));
    h.timebar(ir,2) = line(doy(tplt([1 1])), h.ax(ir,2).YLim, 'parent', h.ax(ir,2), 'color', rgb('gold'));
end

% Legends (have to be added after time bar to prevent additions)

h.leg(1) = legend([h.ts(1).lndaily, h.ts(1).lnsummer, h.ts(1).lnfc], ...
    ["Hindcast (daily)", "Hindcast (Jul 1)", "Forecast"], ...
    'numcolumns', 2, 'location', 'northeast', 'box', 'off');
h.leg(1).Position(2) = h.ax(1,1).Position(2)+h.ax(1,1).Position(4);
h.leg(1).Position(1) = h.ax(1,1).Position(1)+h.ax(1,1).Position(3)-h.leg(1).Position(3);

h.leg(2) = legend(h.doy(1).lndoy(end-3:end), string(2022:2025), ...
    'numcolumns', 2, 'location', 'northwest', 'box', 'off');
h.leg(2).Position(2) = h.ax(1,2).Position(2)+h.ax(1,2).Position(4);
h.leg(2).Position(1) = h.ax(1,2).Position(1);

% Animate

savegif = true;
outfol = '~/Documents/Conferences/202505_SpringPEEC/slides/';

for it = 1:length(tplt)

    if tplt(it) <= datetime(2025,4,11)
        isfc = false;
        MAn = readmom6mapslice(C, tplt(it), ...
        'vars', vplt(iv), 'vartype', 'anomaly');
        MVa = readmom6mapslice(C, tplt(it), ...
            'vars', vplt(iv), 'vartype', 'value');
    else
        isfc = true;
        MFc = readmom6mapslice(C, tplt(it), ...
            'vars', vplt(iv), 'vartype', 'fcpersist');
    end

    h.m(1).pc.CData(1:end-1,1:end-1) = MAn.(vplt{iv});
    if isfc
        h.m(2).pc.CData(1:end-1,1:end-1) = MFc.(vplt{iv});
    else
        h.m(2).pc.CData(1:end-1,1:end-1) = MVa.(vplt{iv});
    end

    if isfc
        h.max(1).Title.String{2} = 'Anomaly: Apr 11, 2025';
        h.max(2).Title.String = sprintf('Prediction: %s', string(tplt(it)));
    else
        h.max(1).Title.String{2} = sprintf('Anomaly: %s', string(tplt(it)));
        h.max(2).Title.String = sprintf('Hindcast: %s', string(tplt(it)));
    end

    for ir = 1:4
        h.timebar(ir,1).XData = tplt([it it]);
        h.timebar(ir,2).XData = doy(tplt([it it]));
    end

    if savegif
        if it == 1
            gif(fullfile(outfol, sprintf('animated_map_index_%s.corrected.gif',vplt)), 'frame', h.fig); %, 'resolution', 150);
        else
            gif;
        end
    else
        drawnow;
    end

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