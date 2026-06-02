%% Anomaly plots and indices slides
%| label: fig-peec25-indices-and-maps
%| fig-cap: From right to left: variable values over time (grey line = daily hindcast values,
%  red line = persistence forecast values, green circles = July 1 values,
%  red circle = July 1 persistence forecast values); variable values versus
%  day of year, highlighting the most recent 4 years with forecast values in the red dashed line; and maps of the (top)
%  spatial anomaly on the last simulated date and (bottom) predicted value
%  based on this anomaly overlaid on July 1 climatological conditions.
%| code-fold: true

% Setup

savetofile = false;
outfol = '~/Documents/Conferences/202605_SpringPEEC/slides/';

% Dataset details

C = cefiportalopts('portalpath', '~/Documents/Data/CEFI/ak_cefiportal/', ...
                   'region', 'nep', ...
                   'release', 'e202604', ...
                   'freq', 'daily', ...
                   'yyyymmdd', '20240101', ...
                   'subdomain', 'iq0-342jq446-743'); 

yr = 2026;

% Regions to plot (survey region masks)

regnum = [98 143 52 47]; % EBS, NBS, AI, GOA
regname = ["EBS", "NBS", "AI", "GOA"];
reglong = ["Southeast Bering Sea", "Northern Bering Sea", "Aleutian Islands", "Gulf of Alaska"];

% Variable setup: which ones do we plot?

vplt = ["tob", "tos", "siconc", "cpool2p0", "cpool0p0"];

% Read surveyavg indices from file

Tbl = readindexdata(C.setopts('subdomain', 'aksvyreg'), ...
    'regnum', regnum, ...
    'vars', vplt);

tfc = Tbl.Fc.t(1);      % forecast start date
tpr = datetime(yr,7,1); % summer prediction date
tfc.Format = 'MMM dd, uuuu';
tpr.Format = 'MMM dd, uuuu';

% Remove leap days from climatology

isleap = month(Tbl.Clim.t) == 2 & day(Tbl.Clim.t) == 29;
Tbl.Clim = Tbl.Clim(~isleap,:);

% Map setup for anomaly/forecast maps

A = akmapprep('cpopts', C, 'maskvar', 'mask_esr_area', 'maskval', [1 3 4], ...
    'latbuffer', 0.01, 'lonbuffer', 0.01);

nreg = length(regnum);
tilesz = [nreg 5];

% Variable info table

V = ak_variable_info;

% Read map slice data

MapAnom = readmom6mapslice(C, Tbl.Fc.t(1), ...
    'vars', setdiff(vplt, ["cpool2p0", "cpool0p0"]), 'vartype', 'anomaly');
MapFc = readmom6mapslice(C, datetime(yr,7,1), ...
    'vars', setdiff(vplt, ["cpool2p0", "cpool0p0"]), 'vartype', 'fcpersist');
MapAnom.cpool2p0 = MapAnom.tob;
MapAnom.cpool0p0 = MapAnom.tob;
MapFc.cpool2p0 = MapFc.tob;
MapFc.cpool0p0 = MapFc.tob;close

% Slide figures

for iv = 1:length(vplt)

    % Figure and tiling setup

    h = struct;
    h.fig = figure('color', 'w');
    h.fig.Position(3:4) = [1300 600]; % slide size
    h.t = tiledlayout(tilesz(1), tilesz(2));
    h.t.TileSpacing = 'loose';
    h.t.Padding = 'compact';

    % Plot index timeseries

    for ir = 1:length(regnum)
        h.ax(ir,1) = nexttile(tilenum(h.t,ir,1), [1 2]);
        h.ts(ir) = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ir, ...
            'style', 'date_vs_val_annual_and_daily');

        h.ax(ir,2) = nexttile(tilenum(h.t,ir,3));
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
    
    h.leg(2) = legend(h.doy(1).lndoy(end-3:end), string((-3:0)+yr), ...
        'numcolumns', 2, 'location', 'northwest', 'box', 'off');
    h.leg(2).Position(2) = h.ax(1,2).Position(2)+h.ax(1,2).Position(4);
    h.leg(2).Position(1) = h.ax(1,2).Position(1);

    % Plot maps

    h.max(1) = nexttile(tilenum(h.t,1,4), [2 2]);
    h.m(1) = plotmom6akmap(A, MapAnom, ...
        'var', vplt{iv}, 'vartype', 'anomaly', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    h.max(2) = nexttile(tilenum(h.t,3,4), [2 2]);
    h.m(2) = plotmom6akmap(A, MapFc, ...
        'var', vplt{iv}, 'vartype', 'value', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});

    % Prettify and label maps

    set([h.max], 'fontsize', 6, 'box', 'on');

    title(h.max(1), {V{vplt{iv}, 'labelname'}, "Anomaly: "+string(tfc)}, 'fontsize', 12);
    title(h.max(2), "Prediction: "+string(tpr), 'fontsize', 12);

    h.m(1).cb.Position(1) = 0.95;
    h.m(2).cb.Position(1) = 0.95;
    set([h.m.cb], 'axislocation', 'out', 'fontsize', 10, 'tickdir', 'out');

    % Save to file

    if savetofile
        export_fig(h.fig, fullfile(outfol, sprintf('peec_persisfc_%s', vplt{iv})), '-png', '-r150', '-nocrop');
    end
end


%% Animated version
%
% Maps of surface and bottom temperature (anomaly and value/prediction)
% with sea ice fraction contours.  Day-of-year index plots of
% surface/bottom provides historical context.

% Time blocks to animate separately (allows me to step through things when
% presenting)

tani = datetime([...
    2025  4 11
    2025  7  1
    2025  9  6
    2025 11  1
    2026  1 16
    2026  2  4
    2026  2 23
    2026  3 24
    2026  4 29
    2026  7  1]);
nani = length(tani)-1;

loop = 0; % loop count for gif files

% Create animated gifs

for iani = 1:nani

    % Times to include

    tplt = tani(iani):tani(iani+1);
    tplt.Format = 'MMM dd, uuuu';
    
    % File setup

    tmp = tplt([1 end]);
    tmp.Format = 'uuuuMMdd';
    outname = fullfile(outfol, sprintf('animated_map_index_tempice_%s-%s.gif', tmp));
    
    if savegif && exist(outname, 'file')
        continue
    end
    
    % Find where these time indices are in the existing hindcast and
    % forecast tables

    idxhc = interp1(Tbl.Hc.t, 1:length(Tbl.Hc.t), tplt, 'nearest');
    idxfc = interp1(Tbl.Fc.t, 1:length(Tbl.Fc.t), tplt, 'nearest');
    
    % Figure and tiling setup
    
    h = struct;
    h.fig = figure('color', 'w');
    h.fig.Position(3:4) = [1300 600]*1.25; % slide size
    h.t = tiledlayout(tilesz(1), tilesz(2));
    h.t.TileSpacing = 'loose';
    h.t.Padding = 'compact';
    
    % Plot index timeseries
    
    for ir = 1:length(regnum)

        % Temperature (surface and bottom)

        h.ax(ir,1) = nexttile(tilenum(h.t,ir,1));
        h.doyb(ir) = plotsurveyregionindex(Tbl, "tob", 'col', ir, ...
            'style', 'doy_vs_val_lines');
        h.doys(ir) = plotsurveyregionindex(Tbl, "tos", 'col', ir, ...
            'style', 'doy_vs_val_lines');
    
        % Add markers to indicate time step of maps

        h.mark(ir,1) = line(doy(Tbl.Hc.t(idxhc(1))),  Tbl.Hc.tob(idxhc(1),ir), 'marker', 'v', 'linestyle', 'none');
        h.mark(ir,2) = line(doy(Tbl.Hc.t(idxhc(1))),  Tbl.Hc.tos(idxhc(1),ir), 'marker', '^', 'linestyle', 'none');
    
        % Sea ice on right axis

        yyaxis(h.ax(ir,1), 'right');
        h.doyi(ir) = plotsurveyregionindex(Tbl, "siconc", 'col', ir, ...
            'style', 'doy_vs_val_lines');
    
        h.mark(ir,3) = line(doy(Tbl.Hc.t(idxhc(1))),  Tbl.Hc.siconc(idxhc(1),ir), 'marker', 'o', 'linestyle', 'none');
    end
    
    % Prettify timeseries
    
    set(h.ax, 'ylim', V{'tos', 'limsvyavg'}, 'fontsize', 8, 'box', 'off');
    set(h.ax(1:end-1), 'xticklabel', '');
    labelaxes(h.ax(:,1), compose("%s (%s)", reglong', regname'), 'northwestoutsideabove', 'fontweight', 'b');
    
    htmp = [h.ax.YAxis]';
    set(htmp(:,2), 'Limits', [0 1], 'color', rgb('navy blue'));
    set(htmp(:,1), 'Limits', V{'tos', 'limsvyavg'});
    htmp(1,1).Label.String = 'Temperature';
    htmp(1,2).Label.String = 'Ice fraction';
    
    set([h.doyi.lndoy], 'marker', 'none', 'linestyle', '-'); % ugh, why does yyaxis *do* this?
    
    set(h.mark, 'markerfacecolor', 'y', 'markeredgecolor', 'k');
    
    set([h.doyb.lndoymean h.doys.lndoymean h.doyi.lndoymean], 'linestyle', '-')
    
    tmp = [h.doyb.lndoy]; set(tmp(1:end-2,:), 'visible', 'off');
    tmp = [h.doys.lndoy]; set(tmp(1:end-2,:), 'visible', 'off');
    tmp = [h.doyi.lndoy]; set(tmp(1:end-2,:), 'visible', 'off');
    
    % Plot maps
    
    MAn = readmom6mapslice(C, tplt(1), ...
            'vars', ["tos", "tob", "siconc"], 'vartype', 'anomaly');
    MVa = readmom6mapslice(C, tplt(1), ...
        'vars', ["tos", "tob", "siconc"], 'vartype', 'value');
    
    % ... surface temp
    
    h.max(1,1) = nexttile(tilenum(h.t,1,2), [2 2]);
    h.m(  1,1) = plotmom6akmap(A, MAn, ...
        'var', 'tos', 'vartype', 'anomaly', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    h.max(2,1) = nexttile(tilenum(h.t,3,2), [2 2]);
    h.m(  2,1) = plotmom6akmap(A, MVa, ...
        'var', 'tos', 'vartype', 'value', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    % ... bottom temp
    
    h.max(1,2) = nexttile(tilenum(h.t,1,4), [2 2]);
    h.m(  1,2) = plotmom6akmap(A, MAn, ...
        'var', 'tob', 'vartype', 'anomaly', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    h.max(2,2) = nexttile(tilenum(h.t,3,4), [2 2]);
    h.m(  2,2) = plotmom6akmap(A, MVa, ...
        'var', 'tob', 'vartype', 'value', ...
        'mapprops', {'mticks', 'b', 'pticks', 'l'});
    
    % ... ice contours
    
    tmpsi = h.m(1).pc.CData;
    tmpsi(1:end-1,1:end-1) = MVa.siconc;
    
    [~,h.c(1,1)] = contour(h.max(1,1), h.m(1).pc.XData, h.m(1).pc.YData, tmpsi, [0.1 0.1], 'color', rgb('navy blue'));
    h.c(1,2) = copyobj(h.c(1,1), h.max(1,2));
    h.c(2,1) = copyobj(h.c(1,1), h.max(2,1));
    h.c(2,2) = copyobj(h.c(1,1), h.max(2,2));
    
    % Prettify and label maps
    
    set([h.max], 'fontsize', 6, 'box', 'on');
    
    % Adjust colormap to highlight cold pool values
    
    cmap = h.max(2,1).Colormap;
    cedge = linspace(h.max(2,1).CLim(1), h.max(2,1).CLim(2), size(cmap,1)+1);
    cmap(cedge < 1.98,:) = brighten(cmap(cedge < 1.98,:), 0.8);
    set(h.max(2,:), 'colormap', cmap);
    
    % Add titles to maps indicating date and variable type
    
    title(h.max(1,1), {V{'tos', 'labelname'}, "Anomaly: "+string(tplt(1))}, 'fontsize', 12);
    title(h.max(2,1), "Value: "+string(tplt(1)), 'fontsize', 12);
    title(h.max(1,2), {V{'tob', 'labelname'}, "Anomaly: "+string(tplt(1))}, 'fontsize', 12);
    title(h.max(2,2), "Value: "+string(tplt(1)), 'fontsize', 12);
    
    % Colorbars for maps

    h.m(1,1).cb.Position(1) = 0.95;
    h.m(2,1).cb.Position(1) = 0.95;
    set([h.m(:,2).cb], 'visible', 'off');
    set([h.m.cb], 'axislocation', 'out', 'fontsize', 10, 'tickdir', 'out');
    
    % Legend for time series

    h.leg = legend(h.doyb(1).lndoy(end-3:end), string((-3:0)+yr), ...
        'numcolumns', 2, 'location', 'northwest', 'box', 'off');
    h.leg.Position(2) = h.ax(1).Position(2)+h.ax(1).Position(4)+.03;
    h.leg.Position(1) = h.ax(1).Position(1);

    % Animate
    
    savegif = true;
    
    for it = 1:length(tplt)

        % Read data for this time slice
    
        if tplt(it) <= tfc
            isfc = false;
            MAn = readmom6mapslice(C, tplt(it), ...
            'vars', ["tos", "tob", "siconc"], 'vartype', 'anomaly');
            MVa = readmom6mapslice(C, tplt(it), ...
                'vars', ["tos", "tob", "siconc"], 'vartype', 'value');
        else
            isfc = true;
            MFc = readmom6mapslice(C, tplt(it), ...
                'vars', ["tos", "tob", "siconc"], 'vartype', 'fcpersist');
        end
    
        % Chane map color data

        h.m(1,1).pc.CData(1:end-1,1:end-1) = MAn.tos;
        h.m(1,2).pc.CData(1:end-1,1:end-1) = MAn.tob;
    
        tmpsi = h.m(1).pc.CData;
        if isfc
            h.m(2,1).pc.CData(1:end-1,1:end-1) = MFc.tos;
            h.m(2,2).pc.CData(1:end-1,1:end-1) = MFc.tob;
            tmpsi(1:end-1,1:end-1) = MFc.siconc;
        else
            h.m(2,1).pc.CData(1:end-1,1:end-1) = MVa.tos;
            h.m(2,2).pc.CData(1:end-1,1:end-1) = MVa.tob;
            tmpsi(1:end-1,1:end-1) = MVa.siconc;
        end
    
        % Recreate ice contours

        delete(h.c);
    
        [~,h.c(1,1)] = contour(h.max(1,1), h.m(1).pc.XData, h.m(1).pc.YData, tmpsi, [0.01 0.01], 'color', rgb('navy blue'));
        h.c(1,2) = copyobj(h.c(1,1,1), h.max(1,2));
        h.c(2,1) = copyobj(h.c(1,1,1), h.max(2,1));
        h.c(2,2) = copyobj(h.c(1,1,1), h.max(2,2));
    
        % Update map titles
    
        if isfc
            for iax = 1:2
                h.max(1,iax).Title.String{2} = sprintf('Anomaly: %s', string(tfc));
                h.max(2,iax).Title.String = sprintf('Prediction: %s', string(tplt(it)));
            end
        else
            for iax = 1:2
                h.max(1,iax).Title.String{2} = sprintf('Anomaly: %s', string(tplt(it)));
                h.max(2,iax).Title.String = sprintf('Hindcast: %s', string(tplt(it)));
            end
        end
    
        % Update timeseries markers

        for ir = 1:4
            if isfc
                h.mark(ir,1).XData = doy(Tbl.Fc.t(idxfc(it)));
                h.mark(ir,1).YData = Tbl.Fc.tos(idxfc(it),ir);
        
                h.mark(ir,2).XData = doy(Tbl.Fc.t(idxfc(it)));
                h.mark(ir,2).YData = Tbl.Fc.tob(idxfc(it),ir);
        
                h.mark(ir,3).XData = doy(Tbl.Fc.t(idxfc(it)));
                h.mark(ir,3).YData = Tbl.Fc.siconc(idxfc(it),ir);
            else
                h.mark(ir,1).XData = doy(Tbl.Hc.t(idxhc(it)));
                h.mark(ir,1).YData = Tbl.Hc.tos(idxhc(it),ir);
        
                h.mark(ir,2).XData = doy(Tbl.Hc.t(idxhc(it)));
                h.mark(ir,2).YData = Tbl.Hc.tob(idxhc(it),ir);
        
                h.mark(ir,3).XData = doy(Tbl.Hc.t(idxhc(it)));
                h.mark(ir,3).YData = Tbl.Hc.siconc(idxhc(it),ir);
            end
        end

        % Save to file (or update screen render)
    
        if savegif
            if it == 1
                gif(outname, 'frame', h.fig, 'LoopCount', loop); %, 'resolution', 150);
            else
                gif;
            end
        else
            drawnow;
        end
    
    end

end




%%
%| label: fig-peec26-cluster-allESR
%| fig-cap: Cluster analysis based on surface and bottom temperature averaged over Apr 1-30 across
%| the three Alaska ESR regions.
%| code-fold: true
%| code-eval: false

[hclus, Cdata] = cluster_on_anomalies(...
    [4 1], ...
    'cpopts', C, ...
    'yr', 1993:2026, ...
    'maskval', [1 3 4], ...
    'ndays', 30, ...
    'verbose', false);

if savetofile
    export_fig(hclus.fig, fullfile(outfol, 'cluster_tostob_0401-4030'), '-png', '-r150');
end


%% TODO: more clusters

% ESR subregions
% AI:  2 7 1
% GOA: 6 5
% EBS: 4
% NBS: 3

for ir = [6 5]

    [hclus, Cdata] = cluster_on_anomalies(...
        [4 1], ...
        'cpopts', C, ...
        'yr', 1993:2026, ...
        'maskvar', "mask_esr_subarea", ...
        'maskval', ir, ...
        'ndays', 30, ...
        'verbose', false, ...
        'Cdata', Cdata);
end

% ESR regions

for ir = [1 3 4]

    [hclus, Cdata] = cluster_on_anomalies(...
        [4 1], ...
        'cpopts', C, ...
        'yr', 1993:2026, ...
        'maskvar', "mask_esr_area", ...
        'maskval', ir, ...
        'ndays', 30, ...
        'verbose', false, ...
        'Cdata', Cdata);
end

% Survey regions

for ir = [47 52 98 143]

    [hclus, Cdata] = cluster_on_anomalies(...
        [4 1], ...
        'cpopts', C, ...
        'yr', 1993:2026, ...
        'maskvar', "mask_survey_area", ...
        'maskval', ir, ...
        'ndays', 30, ...
        'verbose', false, ...
        'Cdata', Cdata, ...
        'vars', "tob");
end

% Whole region

[hclus, Cdata] = cluster_on_anomalies(...
        [4 1], ...
        'cpopts', C, ...
        'yr', 1993:2026, ...
        'maskvar', Grd.wet == 1, ...
        'ndays', 30, ...
        'verbose', false, ...
        'Cdata', Cdata);



