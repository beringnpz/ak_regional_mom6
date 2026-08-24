% Bering10K ROMS + BEST_NPZ vs MOM6-NEP-COBALTv3

%% Map of both domains, with insets showing vertical layers and bathymetry differences

% ROMS Bering10K and MOM6-NEP model grid files (public datasets)

romsgrid = 'https://data.pmel.noaa.gov/aclim/thredds/dodsC/ancillary/Bering10K_extended_grid.nc';
Rgrd = ncstruct(romsgrid);

nepstatic = 'http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northeast_pacific/full_domain/hindcast/daily/raw/r20260701/ocean_static.nc';
Mgrd = ncstruct(nepstatic);

% Map limits: based on larger MOM6 grid

latlim = minmax(Mgrd.geolat_c);
lonlim = minmax(Mgrd.geolon_c);

% Calculate map outline

[xm,ym] = mask2poly(Mgrd.geolon_c, Mgrd.geolat_c, Mgrd.wet==1);
[xr,yr] = mask2poly(Rgrd.lon_psi, Rgrd.lat_psi, Rgrd.mask_rho(2:end-1,2:end-1));

% Color setup

cmap1 = cmocean('-balance');
cmap2 = cptcmap('Paired_12');
cmap2 = cmap2([2 1 6 5],:);
%%
%-------------
% Primary map
%-------------

% Boxed map

h = struct('fig', figure);
h.fig.Position(3:4) = [800 600]*1.5;
marb = 0.02; marl = 0.05;
h.ax = axes('position', [marl marb 1-marl 1-marb]);
h.bm = boxworldmap(latlim, lonlim, 'mticks', 'b', 'pticks', 'l');

padend = @(x) [x nan(size(x,1),1); nan(1,size(x,2)+1)];

% MOM6 grid

h.mask(1) = pcolorm(Mgrd.geolat_c, Mgrd.geolon_c, padend(Mgrd.deptho));

% ROMS grid

rbathy = -Rgrd.h;
rbathy(Rgrd.mask_rho == 0) = NaN;
h.mask(2) = plotromsrho(Rgrd, rbathy, false);

set(gca, 'clim', [-1 1]*7000, 'colormap', cmap1);

% Add outline of ocean grid cell polygons

h.bor(1) = plotm(ym,xm);
h.bor(2) = plotm(yr,xr);
set(h.bor, {'color'}, num2cell(cmap2([1 3],:),2));

set(h.bor, 'linewidth', 1.5);
set(h.mask, 'facealpha', 0.5);

% Add land borders

[blat, blon] = arrayfun(@borders, ["united states", "canada", "russia", "mexico"], 'uni', 0);
[bx,by] = projfwd(getm(h.bm.ax), [blat{:}], [blon{:}]);
plot(polyshape(bx,by),'facecolor', 'k', 'edgecolor', 'none');

% Transect lines across Bering Shelf and then Aleutians

Tsct = struct;
Tsct.latlon = [...
       58.929      -165.09
       55.321      -173.56
       49.538      -169.78
    ];
% Tsct.latlon = [Rgrd.lat_rho(:,89) Rgrd.lon_rho(:,89)];

[Tsct.lat, Tsct.lon] = interpm(Tsct.latlon(:,1), Tsct.latlon(:,2), km2deg(5));
Tsct.dx = distance(Tsct.lat(1:end-1), Tsct.lon(1:end-1), Tsct.lat(2:end), Tsct.lon(2:end), referenceEllipsoid('earth', 'km'));
Tsct.x = cumsum([0; Tsct.dx]);
[d,id] = pdist2([Tsct.lat Tsct.lon], Tsct.latlon(2,:), 'euclidean', 'smallest', 1);

mask = Mgrd.wet == 1;
fs = scatteredInterpolant(Mgrd.geolon(mask), Mgrd.geolat(mask), Mgrd.deptho(mask), 'linear');
Tsct.zbot(:,1) = fs(wrapTo360(Tsct.lon), Tsct.lat); 
fs = scatteredInterpolant(Rgrd.lon_rho(:), Rgrd.lat_rho(:), Rgrd.h(:), 'linear');
Tsct.zbot(:,2) = fs(wrapTo360(Tsct.lon), Tsct.lat); 
mask = ~isnan(Mgrd.geolon);

plotm(Tsct.lat, Tsct.lon, '--k');

%----------------
% Vertical layers
%----------------

zmin = 0.03;

% MOM6 layers (as polygons)

xwater = [Tsct.x(1); Tsct.x(end); Tsct.x(end:-1:1,1)];
ywater = [zmin        ; zmin          ; Tsct.zbot(end:-1:1,1)];
pwater = polyshape(xwater, ywater);

xlayer = Tsct.x([1 1 end end 1]);
ylayer = [Mgrd.zi(1:end-1) Mgrd.zi(2:end) Mgrd.zi(2:end) Mgrd.zi(1:end-1) Mgrd.zi(1:end-1)];

for ii = size(ylayer, 1):-1:1
    tmp = polyshape(xlayer', ylayer(ii,:));
    pm(ii) = tmp.intersect(pwater);
end

% ROMS layers (as polygons)

[~,zw] = calcromsz(Tsct.zbot(:,2), -zmin, 30);
zw = permute(zw, [1 3 2]);

for ii = size(zw,2)-1:-1:1
    pr(ii) = polyshape([Tsct.x; Tsct.x(end:-1:1)], -[zw(:,ii); zw(end:-1:1,ii+1)]);
end

% Create new axis

drawnow;
pos = plotboxpos(h.bm.ax);
h.axt = axes('position', [pos(1)+pos(3)*0.01 pos(2)+pos(4)*0.05 pos(3)*0.4, pos(4)*0.35]);

% Plot layers

h.bottom = plot(h.axt, Tsct.x, Tsct.zbot);
set(h.bottom, {'color'}, num2cell(cmap2([1 3],:),2));

hold(h.axt, 'on');
h.lay{1} = plot(h.axt, pm, 'facecolor', cmap2(1,:), 'edgecolor', 'none');
h.lay{2} = plot(h.axt, pr, 'facecolor', cmap2(3,:), 'edgecolor', 'none');
set(h.lay{1}(2:2:end), 'facecolor', 'w');
set(h.lay{2}(end-1:-2:1), 'facecolor', 'w');
    
% Configure axis

set(h.axt, 'ydir', 'reverse', 'xdir', 'reverse', ...
   'yaxisloc', 'right', 'fontsize', 8, 'xlim', minmax(Tsct.x), ...
   'xaxisloc', 'top', 'yscale', 'log', 'ylim', [.2 7000]);
gridxy(Tsct.x(id), [], 'color', 'k', 'linestyle', ':', 'parent', h.axt);

xlabel(h.axt, 'Distance along dashed transect (km)');
ylabel('Depth (m)');

%----------------------------------
% Bathymetry difference, MOM6-ROMS
%----------------------------------

% Regrid MOM6 values onto ROMS grid for comparison

fs = scatteredInterpolant(Mgrd.geolon(mask), Mgrd.geolat(mask), Mgrd.deptho(mask), 'linear');
zmom6onroms = fs(Rgrd.lon_rho, Rgrd.lat_rho);
dbathy = zmom6onroms - Rgrd.h;

% Add axis in upper right

h.axd = axes('position', [pos(1)+pos(3)-pos(3)*0.4 pos(2)+pos(4)-pos(4)*0.4 pos(3)*0.4, pos(4)*0.4]);

% Plot map

plotromsrho(Rgrd, dbathy);
cmocean('balance', 'pivot', 0);
setm(h.axd, 'ffacecolor', 'w', 'fontsize', 7, 'flinewidth', 1);
set(h.axd, 'clim', [-1000 1000], 'colormap', cmocean('diff'));

% Colorbar

h.cb = colorbar('north');
h.cb.Position = [h.cb.Position(1)+h.cb.Position(3)*0.25 ...
                 h.cb.Position(2)-0.07 ...
                 h.cb.Position(3)*0.5 ...
                 h.cb.Position(4)];
xlabel(h.cb, "Modeled bathymetry difference (m)");
set(h.cb, 'axislocation', 'out');
h.an(1) = annotation('textbox', h.cb.Position, 'string', 'ROMS deeper', 'horiz', 'left'); %, 'vert','base','horiz','left')
h.an(2) = annotation('textbox', h.cb.Position, 'string', 'MOM6 deeper', 'horiz', 'right'); %, 'vert','base','horiz','left')
set(h.an, 'color', rgb('light gray'), 'fontsize', 8, 'edgecolor', 'none')

%% Polygon-based skill analysis: setup

%% ... default polygon setup

[Geo, Sgrd, Tri] = ak_polygons;
reg = fieldnames(Geo);
nreg = length(reg);

svyid = cellfun(@(x) unique([Geo.(x).SURVEY_DEFINITION_ID]), reg);

tmp = cell2mat(struct2cell(structfun(@(x) x.verts, Tri, 'uni', 0)));
latlim = minmax(tmp(:,2));
lonlim = minmax(wrapTo360(tmp(:,1)));

%% ... default map

worldmap(latlim, lonlim);
for ir = 1:nr
    plotpolyfv(Tri.(reg{ir}));
end

bordersm('counties', 'facecolor', ones(1,3)*0.8, 'edgecolor', 'none');

%% ... how does the area of each polygon compare?

worldmap(latlim, lonlim);
for ir = 1:nr
    plotpolyfv(Tri.(reg{ir}), 'val', [Geo.(reg{ir}).AREA_M2]/1e6);
end

bordersm('counties', 'facecolor', ones(1,3)*0.8, 'edgecolor', 'none');

%%  Survey and survey-replicated data

% Read data from most recent AKFIN query

Tmp = readtable('svy_temp.csv');
Tmp.LONGITUDE = wrapTo360(Tmp.LONGITUDE);
Tmp = Tmp(~(Tmp.STATION == ""),:);
% Tmp.STATION{Tmp.STATION == ""} = 'null';

% Split by region

for ii = 1:nr
    Svy.(reg{ii}) = Tmp(Tmp.SURVEY_DEFINITION_ID == svyid(ii),:);
end

% Dataset info

reginfo = table(svyid, reg, 'variablenames', ["id","name"]);

vname = ["bt", "sst"];
dataname = [...
    "GEAR_TEMPERATURE_C", "SURFACE_TEMPERATURE_C"
    "mom6nep_r20250912_tob", "mom6nep_r20250912_tos"
    "romsb10k_K20_CORECFS_temp_bottom5m", "romsb10k_K20_CORECFS_temp_surface5m"
    "glorys_q20260716_bottomT", "glorys_q20260716_thetao"];
nr = height(reginfo);
nv = length(vname);
nd = size(dataname,1);



%% ... Collect survey[rep] data into polygon x year x dataset x variable arrays

for ii = 1:nr
    x = reginfo.name{ii};

    % Indetify points by polygon and year
    [yr, ~,iyr] = unique(Svy.(x).YEAR);
    if ismember(x, {'ebs','nbs'})
        [tf,ip] = ismember(Svy.(x).STATION, [Geo.(x).STATION]);
        if ~all(tf)
            error('Mismatch in ebs/nbs?');
        end
    else
        [tf,ip] = ismember(Svy.(x).STRATUM, [Geo.(x).AREA_ID]);

        if ~all(tf)
            [xtmp,ytmp] = projfwd(Geo.(x)(1).Shape.ProjectedCRS, Svy.(x).LATITUDE(~tf), Svy.(x).LONGITUDE(~tf));
            
            [unq,~,iunq] = unique(Svy.(x).STRATUM(~tf));
            pidx = nan(nnz(~tf),1);
            for iu = 1:length(unq)
                ispolysub = floor([Geo.(x).AREA_ID]) == unq(iu);
                ind = find(ispolysub);
                [in,index] = inpolygons(xtmp(iunq==iu), ytmp(iunq==iu), [Geo.(x)(ispolysub).X], [Geo.(x)(ispolysub).Y]);
                pidx(iunq==iu) = ind(cell2mat(index));
            end

            ip(~tf) = pidx;

        end

    end

    np = length(Geo.(x));
    ny = length(yr);

    [gid, ~, g] = unique([ip iyr], 'rows');
    idx = sub2ind([np ny], gid(:,1), gid(:,2));

    DataP.(x).val = cell(np,ny,nd,nv);
    DataP.(x).year = yr;

    for id = 1:nd
        for iv = 1:nv
            tmp = cell(np,ny);
            [tmp{:}] = deal([NaN]);
            datagrp = splitapply(@(x) {x}, Svy.(x).(dataname{id,iv}), g);
            tmp(idx) = datagrp;

            DataP.(x).val(:,:,id,iv) = tmp;

        end
    end

end