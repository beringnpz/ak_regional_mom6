function cluster_on_anomalies(simname, mmdd, opt)

arguments
    simname {mustBeTextScalar}
    mmdd (1,2) {mustBeInteger}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
    opt.yr {mustBeInteger, mustBeVector} =(1993:year(datetime('today')))
    opt.vars {mustBeText} =["tob", "tos"]
    opt.maskvar ="mask_esr_area"
    opt.maskval =[1 2 3 4];
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

for iy = nyr:-1:1
    fprintf('Reading data: %d\n', opt.yr(iy));

    % Check for files for this year

    fglob = fullfile(opt.datafol, simname, "Level3", simname+"_daily_anomaly_" + opt.yr(iy) + "*.nc");
    fanom = dir(fglob);
    nfile = length(fanom);
    if nfile < 1
        error('No anomaly files found for %d; adjust yr input accordingly', opt.yr(iy));
    end
    fanom = fullfile({fname.folder}, {fname.name});

    % Find time index corresponding to clustering target date

    t = ncdateread(fanom, 'time');
    tidx = interp1(t, 1:length(t), datetime(opt.yr(iy), mmdd(1), mmdd(2)), 'nearest');

    Cdata(iy) = ncstruct(fanom, struct('time', [tidx 1 1]), opt.vars{:});

end

return

% TODO, old code below for reference

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