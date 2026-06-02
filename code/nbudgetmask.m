%% Mask-building: classify entire domain

rebuildmask = true;
plotmask = false;

Ctest = cefiportalopts('region', 'nep', 'release', 'r20250912', 'freq', 'monthly');
Grd = readcefigridvars(Ctest, ["deptho", "geolat", "geolon", "geolat_c", "geolon_c", "wet"], 'expandname', false);

if rebuildmask

    % Main Bering polygon from marine areas

    Marea = shapeprjread('/Volumes//LaCie2023/NaturalEarth/ne_10m_geography_marine_polys/ne_10m_geography_marine_polys.shp');
    beringareas = {'Bering Sea', 'Gulf of Anadyr''', 'Karaginskiy Gulf', ...
                   'Norton Sound', 'Baird Inlet', 'Bristol Bay'};
    isin = ismember({Marea.name}, beringareas);   
    Marea = Marea(isin);

    [~,mloc] = ismember(beringareas, {Marea.name});

    p = polyshape(wrapTo360([Marea.Lon]), [Marea.Lat]);

    isber = isinterior(p, Grd.geolon(:), Grd.geolat(:)) & Grd.wet(:)==1;
    isber = reshape(isber, size(Grd.deptho));

    % Shelf from model 200-m contour

    % TODO: revisit, doesn't work b/c of NaNs in coordinates

    htmp = figure;
    [cc,hc] = contour(Grd.geolon, Grd.geolat, Grd.deptho, [200 200], 'k');
    C = contourcs(cc, 'cmat');
    [~,imax] = max([C.Length]);
    C = C(imax);
    close(htmp);

    ltshelf = [C.Y Grd.geolat(end,:) Grd.geolat(end:-1:1,end)'];
    lnshelf = wrapTo360([C.X Grd.geolon(end,:) Grd.geolon(end:-1:1,end)']);

    isshelf = inpolygon(Grd.geolon, Grd.geolat, lnshelf, ltshelf) & ...
              isber & Grd.wet == 1 & ...
              Grd.geolat > 54.83; % cut off bit south of unimak

    % Trim shelf between west and east

    iscw = ispolycw(Marea(mloc(1)).Lon, Marea(mloc(1)).Lat);
    [bslat, bslon] = polysplit(Marea(2).Lat, Marea(2).Lon);
    iseast = cellfun(@(x) all(x>0), bslon);
    lteast = bslat{iscw' & iseast};
    lneast = bslon{iscw' & iseast};

    isebs = isshelf & ...
            ~inpolygon(Grd.geolon, Grd.geolat, lneast, lteast) & ...
            Grd.geolon > 170;

    iswbs = isshelf & ~isebs;

    % Specific inlets

    isan = strcmp({Marea.name}, 'Gulf of Anadyr''');
    panadyr = polyshape(wrapTo360([Marea(isan).Lon]), [Marea(isan).Lat]);
    isanadyr = panadyr.isinterior(Grd.geolon(:), Grd.geolat(:)) & Grd.wet(:)==1;
    isanadyr = reshape(isanadyr, size(Grd.deptho));
    isnorton = inpolygon(Grd.geolon, Grd.geolat, ...
        wrapTo360(Marea(mloc(4)).Lon), Marea(mloc(4)).Lat) & Grd.wet==1;

    % Now classify

    masks = {...
        'SEBSinner'     isebs & Grd.geolat < 60.83 & Grd.deptho <= 50
        'SEBSmiddle'    isebs & Grd.geolat < 60.83 & Grd.deptho <= 100 & Grd.deptho > 50
        'SEBSouter'     isebs & Grd.geolat < 60.83 & Grd.deptho > 100
        'NEBSinner'     isebs & ~isanadyr & ~isnorton & Grd.geolat > 60.83 & Grd.deptho <= 50
        'NEBSmiddle'    isebs & ~isanadyr & ~isnorton & Grd.geolat > 60.83 & Grd.deptho <= 100 & Grd.deptho > 50 
        'NEBSouter'     isebs & ~isanadyr & ~isnorton & Grd.geolat > 60.83 & Grd.deptho > 100
        'NortonSound'   isnorton
        'GulfOfAnadyr'  isanadyr
        'WBS'           iswbs
        'Basin'         isber & ~isshelf
        'Chukchi'       Grd.geolat > 65 & ~isebs & Grd.wet == 1};

    masks{2,2} = imfill(masks{2,2}, 'holes') & Grd.wet == 1; % reassign St Paul/George/Matt from inner to middle
    masks{1,2} = masks{1,2} & ~masks{2,2};
%     masks{5,2} = imfill(masks{5,2}, 'holes') & Grd.wet == 1; % same for St. Matt
%     masks{4,2} = masks{4,2} & ~masks{5,2};
    masks{4,2} = imfill(masks{4,2}, 'holes'); % keep deeper bit near Bering Strait in middle
    masks{5,2} = masks{5,2} & ~masks{4,2};
    masks{10,2} = imfill(masks{10,2}, 'holes') & Grd.wet == 1;

    masks = [masks; {'GulfOfAlaska', Grd.wet == 1 & ~any(cat(3, masks{:,2}),3)}];
    
    save('nbudget_masks', 'masks');
else
    masks = load('nbudget_masks.mat');
    masks = masks.masks;
end

nmask = size(masks,1);