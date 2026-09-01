

%% Simulation details

% addpath('../code');

C = cefiportalopts('portalpath', '~/Documents/Data/CEFI/ak_cefiportal', ...
                   'release', 'e202604', ...
                   'subdomain', [0 342 446 743], ...
                   'yyyymmdd', '20240101');

yr = 2026;
esrfol = '~/Documents/Manuscripts/202608_ESR/';

%% Map of bottom temperature, past 20 years

h = plot_bering_bottom_temp_map_5x4(yr, C);

archivepng(h.fig, esrfol, sprintf('btemp_map_%d-%d', [yr-20+1 yr]), ...
    '-r150', '-nocrop');

%% Region-based analysis setup

regnum = [98 143 52 47]; % EBS, NBS, AI, GOA
regname = ["EBS", "NBS", "AI", "GOA"];

% Variable setup: which ones do we plot?

vplt = ["tob", "tos", "cpool2p0", "cpool0p0"];

% Read surveyavg indices from file

Tbl = readindexdata(C.setopts('subdomain', 'aksvyreg'), ...
    'regnum', regnum, ...
    'vars', vplt);

% Remove leap days from climatology

isleap = month(Tbl.Clim.t) == 2 & day(Tbl.Clim.t) == 29;
Tbl.Clim = Tbl.Clim(~isleap,:);

% Extract July 1 values from each year

Jul1 = structfun(@(x) retime(x, datetime(1993:yr,7,1), 'nearest', 'endvalues', NaN), Tbl, 'uni', 0);

% Plot timeseries: ebs

V = ak_variable_info;

%% ... EBS timeseries

ii = 1;
for iv = 1:length(vplt)
        
    h = plotgrid('size', [1 1], 'mar', 0.08, 'mr', 0.01);
    h.fig.Position(3:4) = [6 4]*72;
    h.fig.Color = 'w';

    h.gobj = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ii, 'year', yr);

    title(h.ax, sprintf('%s: %s', regname(ii), V{vplt{iv},'labelname'}));

    archivepng(h.fig, esrfol, sprintf('timeseries_%s_%s',regname{ii}, vplt{iv}), ...
        '-r150', '-nocrop')

end

%% ... NBS timeseries

ii = 2;
for iv = 1:length(vplt)
        
    h = plotgrid('size', [1 1], 'mar', 0.08, 'mr', 0.01);
    h.fig.Position(3:4) = [6 4]*72;
    h.fig.Color = 'w';

    h.gobj = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ii, 'year', yr);

    title(h.ax, sprintf('%s: %s', regname(ii), V{vplt{iv},'labelname'}));

    archivepng(h.fig, esrfol, sprintf('timeseries_%s_%s',regname{ii}, vplt{iv}), ...
        '-r150', '-nocrop')

end

%% ... Bering Sea cluster

% Summer (mean, Jun1-Aug13)

ndays = days(datetime(2026,8,13)-datetime(2026,6,1));

[hclus(1,1), Cdata] = cluster_on_anomalies(...
    [6 1], ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 3, ...
    'ndays', ndays, ...
    'verbose', false);

archivepng(hclus(1,1).fig, esrfol, 'cluster_EBS_summer', '-r300');

% winter-to-summer evolution (mean first week of every month, feb-jul)

[hclus(1,2), Cdata2] = cluster_on_anomalies(...
    [2:7; ones(1,6)]', ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 3, ...
    'ndays', 7, ...
    'verbose', false);

archivepng(hclus(1,2).fig, esrfol, 'cluster_EBS_feb-jul', '-r300');


%% ... Aleutians timeseries

ii = 3;
for iv = 1:2 %length(vplt)
        
    h = plotgrid('size', [1 1], 'mar', 0.08, 'mr', 0.01);
    h.fig.Position(3:4) = [6 4]*72;
    h.fig.Color = 'w';

    h.gobj = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ii, 'year', yr);

    title(h.ax, sprintf('%s: %s', regname(ii), V{vplt{iv},'labelname'}));

    archivepng(h.fig, esrfol, sprintf('timeseries_%s_%s',regname{ii}, vplt{iv}), ...
        '-r150', '-nocrop');

end

%% ... Aleutians cluster

hclus(2,1) = cluster_on_anomalies(...
    [6 1], ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 1, ...
    'ndays', ndays, ...
    'verbose', false, ...
    'Cdata', Cdata);

archivepng(hclus(2,1).fig, esrfol, 'cluster_AI_summer', '-r300');

% winter-to-summer evolution (mean first week of every month, feb-jul)

hclus(2,2) = cluster_on_anomalies(...
    [2:7; ones(1,6)]', ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 1, ...
    'ndays', 7, ...
    'verbose', false, ...
    'Cdata', Cdata2);

archivepng(hclus(2,2).fig, esrfol, 'cluster_AI_feb-jul', '-r300');

%% ... GOA timeseries

ii = 4;
for iv = 1:2 %length(vplt)
        
    h = plotgrid('size', [1 1], 'mar', 0.08, 'mr', 0.01);
    h.fig.Position(3:4) = [6 4]*72;
    h.fig.Color = 'w';

    h.gobj = plotsurveyregionindex(Tbl, vplt{iv}, 'col', ii, 'year', yr);

    title(h.ax, sprintf('%s: %s', regname(ii), V{vplt{iv},'labelname'}));

    archivepng(h.fig, esrfol, sprintf('timeseries_%s_%s',regname{ii}, vplt{iv}), ...
        '-r150', '-nocrop');

end

%% ... GOA cluster

hclus(3,1) = cluster_on_anomalies(...
    [6 1], ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 4, ...
    'ndays', ndays, ...
    'verbose', false, ...
    'Cdata', Cdata);

archivepng(hclus(3,1).fig, esrfol, 'cluster_GOA_summer', '-r300');

% winter-to-summer evolution (mean first week of every month, feb-jul)

hclus(3,2) = cluster_on_anomalies(...
    [2:7; ones(1,6)]', ...
    'cpopts', C, ...
    'yr', 1993:yr, ...
    'maskval', 4, ...
    'ndays', 7, ...
    'verbose', false, ...
    'Cdata', Cdata2);

archivepng(hclus(3,2).fig, esrfol, 'cluster_GOA_feb-jul', '-r300');

%% All regions: spatial evolution of values/anomalies, Jan-Aug

[A, Grd] = akmapprep('cpopts', C, 'maskval', [1 3 4]);

% Data to plot

ttarget = datetime(2026,1:8,1);

pltvars = ["siconc", "tob", "tos",     "siconc",  "tob",     "tos"];
plttype = ["value",  "value", "value", "anomaly", "anomaly", "anomaly"];
pltlbl =  ["Sea ice", "Bottom temp. (\circC)", "SST (\circC)", ...
           "Anomaly: Sea ice", "Anomaly: Bottom temp. (\circC)", "Anomaly: SST (\circC)"];

% Map and axis setup

mapprops = {'mtick', '', 'pticks', ''};

nt = length(ttarget);
np = length(pltvars);

h = plotgrid('size', [nt np], 'sp', 0, 'mar', 0, 'mt', 0.08);
h.fig.Position(3:4) = [800 600]*1.5;
h.fig.Color = 'w';

% Plot

for it = 1:nt
    for iv = 1:np
        if it == 1
            cbloc = 'north';
        else
            cbloc = [];
        end

        Tmp = readmom6mapslice(C, ttarget(it), 'vars', pltvars{iv}, 'vartype', plttype{iv});
        
        axes(h.ax(it,iv));
        h.map{it,iv} = plotmom6akmap(A, Tmp, ...
            'var', pltvars{iv}, ...
            'vartype', plttype{iv}, ...
            'mapprops', mapprops, ...
            'cbloc', cbloc);

        if it == 1
            h.map{it,iv}.cb.Position(2) = 0.94;
            xlabel(h.map{it,iv}.cb, pltlbl{iv});
            set(h.map{it,iv}.cb, 'axisloc', 'out', 'tickdir', 'out');
        end
    end
end

% Label dates

labelaxes(h.ax(:,1), string(ttarget), 'northwest');

% Save to file

archivepng(h.fig, esrfol, sprintf('val_anom_%d', yr), '-r150', '-nocrop');

%% How well did the persistence signal forecast summer values?

[A, Grd] = akmapprep('cpopts', C, 'maskval', [1 3 4]);

ttarget = datetime(2026,6:8,1, 'format', ', yyyy/MM/dd');
nt = length(ttarget);
pltvars = ["tob", "tos"];
pltvname = ["Bottom temp. (\circC)", "SST (\circC)"];
np = length(pltvars);

vartype = ["value", "fcpersis", "anomaly"];
nv = length(vartype);

mapprops = {'mtick', '', 'pticks', ''};

% Plot realized value, predicted value, realized anomaly, and anomaly
% mismatch

h = plotgrid('size', [nt+1, np], 'mar', 0.05, 'sv', 0.05, 'sh', 0.05);
h.fig.Position(3:4) = [900 1062];
h.ax = subgridaxes(h.ax, 2,2);

cbloc = 'east';

for it = 1:nt
    for iv = 1:np

        axes(h.ax(1,1,it+1,iv));
        Tmp1 = readmom6mapslice(C, ttarget(it), 'vars', pltvars{iv}, 'vartype', 'value');
        htmp = plotmom6akmap(A, Tmp1, ...
            'var', pltvars{iv}, ...
            'vartype', 'value', ...
            'mapprops', mapprops, ...
            'cbloc', cbloc);
        reposition(htmp.cb, [0.5 0.8], [0 0], 'e');

        axes(h.ax(1,2,it+1,iv));
        Tmp2 = readmom6mapslice(C, ttarget(it), 'vars', pltvars{iv}, 'vartype', 'fcpersist');
        htmp = plotmom6akmap(A, Tmp2, ...
            'var', pltvars{iv}, ...
            'vartype', 'value', ...
            'mapprops', mapprops, ...
            'cbloc', cbloc);
        reposition(htmp.cb, [0.5 0.8], [0 0], 'e');

        axes(h.ax(2,1,it+1,iv));
        Tmp3 = readmom6mapslice(C, ttarget(it), 'vars', pltvars{iv}, 'vartype', 'anomaly');
        htmp = plotmom6akmap(A, Tmp3, ...
            'var', pltvars{iv}, ...
            'vartype', 'anomaly', ...
            'mapprops', mapprops, ...
            'cbloc', cbloc);
        reposition(htmp.cb, [0.5 0.8], [0 0], 'e');

        axes(h.ax(2,2,it+1,iv));
          
        Tmp4.(pltvars{iv}) = Tmp1.(pltvars{iv}) - Tmp2.(pltvars{iv});
        htmp = plotmom6akmap(A, Tmp4, ...
                    'var', pltvars{iv}, ...
                    'vartype', 'anomaly', ...
                    'mapprops', mapprops, ...
                    'cbloc', cbloc);
        reposition(htmp.cb, [0.5 0.8], [0 0], 'e');

        labelaxes(h.ax(:,:,it+1,iv), ["a)", "c)", "b)", "d)"], 'northwest');

        % Add predicted anomaly

        if it == 1
            fcanom = Tmp2.(pltvars{iv}) - (Tmp1.(pltvars{iv}) - Tmp3.(pltvars{iv}));
            axes(h.ax(2,1,1,iv));
            htmp = plotmom6akmap(A, struct(pltvars{iv},fcanom), ...
                    'var', pltvars{iv}, ...
                    'vartype', 'anomaly', ...
                    'mapprops', mapprops, ...
                    'cbloc', cbloc);
            reposition(htmp.cb, [0.5 0.8], [0 0], 'e');
        end


    end
end

lbls = pltvname + string(ttarget');
labelaxes(h.ax(1,1,2:end,:), lbls, 'northwestoutsideabove');

lblsfc = pltvname + ": forecast anomaly";
labelaxes(h.ax(2,1,1,:), lblsfc, 'northwestoutsideabove')


set(h.ax(2,2,:,:), 'colormap', cmocean('curl'));
delete(h.ax(1,:,1,:));
delete(h.ax(2,2,1,:));

archivepng(h.fig, esrfol, sprintf('fcpersis_vs_hindcast_%d', yr), '-r150');

%% How well did this year's hindcast match the groundfish surveys?

%--------------------
% Setup
%--------------------

% Survey and survey grid polygons

[Geo, Sgrd, Tri] = ak_polygons;

reg = fieldnames(Geo);
nreg = length(reg);

svyid = cellfun(@(x) unique([Geo.(x).SURVEY_DEFINITION_ID]), reg);

% Read data from most recent AKFIN query

% TODO: need to requery AKFIN for 2026 data, change year and release in
% code below

% svycsv = "~/Documents/Manuscripts/2026_MOM6_ak_skill/data/gfsvy_q20260826.csv";
svycsv = "~/Documents/Manuscripts/2026_MOM6_ak_skill/code/svy_temp.csv";

yrsvy = 2024; % TODO: placeholder
Tmp = readtable(svycsv);
Tmp.LONGITUDE = wrapTo360(Tmp.LONGITUDE);
mask = ~(Tmp.STATION == "") & Tmp.YEAR == yrsvy;
Tmp = Tmp(mask,:);

% Split survey data by region

Svy = struct;
for ii = 1:nreg
    Svy.(reg{ii}) = Tmp(Tmp.SURVEY_DEFINITION_ID == svyid(ii),:);
end

% Map setup for this year's survey regions

hassurvey = cellfun(@(x) ~isempty(Svy.(x)), reg);

A = akmapprep('cpopts', C, 'maskvar', 'mask_survey_area', 'maskval', svyid(hassurvey));

% Variables: surface and bottom temp
% Datsets: survey and survey-resampled-MOM6

vcomp = ["GEAR_TEMPERATURE_C",   "mom6nep_r20250912_tob"; ...
         "SURFACE_TEMPERATURE_C" "mom6nep_r20250912_tos"];
[nv,nd] = size(vcomp);

% Rearrange data to polygon x variable x dataset arrays

svydata = struct;
for ir = 1:nreg
    if hassurvey(ir)
        rr = reg{ir};

        svydata.(rr) = nan(length(Geo.(rr)), 2, 2);

        switch rr
            case {'ebs', 'nbs'} % max 1 sample per station, simple lookup
                [tf,loc] = ismember([Geo.(rr).STATION], Svy.(rr).STATION);

                for iv = 1:nv
                    for id = 1:nd
                        svydata.(rr)(tf,iv,id) = Svy.(rr).(vcomp{iv,id})(loc(tf)); 
                    end
                end

            case {'ai', 'goa'} % several stations per area, calculate means

                [~,loc] = ismember(Svy.(rr).STATION, [Sgrd.(rr).STATION]);
                Svy.(rr).AREA_ID = [Sgrd.(rr)(loc).AREA_ID]';
                
                [iunq, unqid] = findgroups(Svy.(rr).AREA_ID);
                nper = splitapply(@length, iunq, iunq);
                
                [tf,loc] = ismember(unqid, [Geo.(rr).AREA_ID]);
                
                for iv = 1:nv
                    for id = 1:nd
                        svydata.(rr)(loc,iv,id) = splitapply(@mean, Svy.(rr).(vcomp{id,iv}), iunq);
                    end
                end
        end
    end
end

%% Plot

h = plotgrid('size', [nd+1 nv], 'sp', 0, 'mar', 0.02);
h.fig.Position(3:4) = [850 1100];

for iv = 1:nv
    for id = 1:nd
        axes(h.ax(id,iv));

        worldmap(A.latlim, A.lonlim);
        for ir = 1:nreg
            if hassurvey(ir)
                h.p(iv,id,ir) = plotpolyfv(Tri.(reg{ir}), 'val', svydata.(reg{ir})(:,iv,id), 'nancol', ones(1,3)*0.8);
                h.s(iv,id,ir) = scatterm(Svy.ebs.LATITUDE, Svy.ebs.LONGITUDE, 2, Svy.ebs.(vcomp{iv,id}), 'filled', 'markeredgecolor', 'k');
            end
        end

    end

    axes(h.ax(nd+1,iv));

    worldmap(A.latlim, A.lonlim);
    for ir = 1:nreg
        if hassurvey(ir)
            h.p(iv,nd+1,ir) = plotpolyfv(Tri.(reg{ir}), 'val', svydata.(reg{ir})(:,iv,2) - svydata.(reg{ir})(:,iv,1), 'nancol', ones(1,3)*0.8);
            h.s(iv,nd+1,ir) = scatterm(Svy.(reg{ir}).LATITUDE, Svy.(reg{ir}).LONGITUDE, 2, Svy.(reg{ir}).(vcomp{iv,2})-Svy.(reg{ir}).(vcomp{iv,1}), 'filled', 'markeredgecolor', 'k');
        end
    end

end

set(h.ax(1:nd,:), 'clim', [-2 15], 'colormap', cmocean('-dense', 8*4));
set(h.ax(nd+1,:), 'clim', [-3 3], 'colormap', cmocean('balance'));


V = ak_variable_info;
cmap = V{'tob','cmap'}{1};
clim = V{'tob','limmap'};
cedge = linspace(clim(1), clim(2), size(cmap,1)+1);
cmap(cedge < 1.98,:) = brighten(cmap(cedge < 1.98,:), 0.8);

set(h.ax(1:nd,:), 'clim', clim, 'colormap', cmap);
set(h.ax(nd+1,:), 'clim', [-3 3], 'colormap', cmocean('balance'));


%%

function reposition(hobj, scale, shift, anchor)
%REPOSITION Scale and/or shift a graphics object's position
%
% Input variables:
%
%   hobj:   handle to graphics object
%   scale:  1x2 array, scale factor to apply to horizontal and vertical
%           positions, respectively
%   shift:  1x2 array: shift to apply to anchor point position, in units
%           used by hobj Position attribute
%   anchor: anchor point indicating corner or side to anchor in place
%           during rescaling  

    arguments
        hobj (1,1)
        scale (1,2) {mustBeNumeric}
        shift (1,2) {mustBeNumeric}
        anchor {mustBeTextScalar, mustBeMember(anchor, {'nw','n','ne','e','se','s','sw','w'})}
    end

    if ~isgraphics(hobj)
        error('First input (hobj) must be valid graphics object');
    end

    pos = get(hobj, 'Position');
    
    w = pos(3)*scale(1);
    h = pos(4)*scale(2);
    
    A = table(...
        [...
            pos(1)          pos(2)+pos(4)
            pos(1)+pos(3)/2 pos(2)+pos(4)
            pos(1)+pos(3)   pos(2)+pos(4)
            pos(1)+pos(3)   pos(2)+pos(4)/2
            pos(1)+pos(3)   pos(2)
            pos(1)+pos(3)/2 pos(2)
            pos(1)          pos(2)
            pos(1)          pos(2)+pos(4)/2], ...
        [...
            0       -h
            -w/2    -h
            -w      -h
            -w      -h/2
            -w      0
            -w/2    0
            0       0
            0       -h/2], ...
        'rownames', {'nw','n','ne','e','se','s','sw','w'}, ...
        'variablenames', {'refxy', 'shift'});
    
    corner = A{anchor,'refxy'} + shift + A{anchor,'shift'};
    newpos = [corner w h];
    
    set(hobj, 'position', newpos);

end