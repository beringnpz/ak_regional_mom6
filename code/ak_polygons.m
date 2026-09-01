function [Geo, Sgrd, Tri] = ak_polygons(opt)
%AK_POLYGONS Create polygons for Alaska strata/survey-based analysis
%
% This function returns structures to support analysis and plotting of the
% Alaska survey strata and survey grid.  It's primarily used for MOM6/ROMS
% comparison to survey datasets.
%
% Output variables:
%
%   Geo:    structure of geographic data structures with polygon data from
%           the survey strata (EBS and NBS) and survey grid (GOA, AI, BSS)
%           datasets.
%
%   Sgrd:   survey grid polygons for all regions.  For EBS and NBS, this is
%           identical to Geo.
%
%   Tri:    face/vertex data for the Geo polygons, used for plotting
%           polygons with scaled color data (see plotpolyfv).  Fields are
%           as follows:
%
%           faces:  nf x 1, face indices
%           verts:  nv x 2, vertex coordinates ([longitude latitude])
%           fvind:  nf x 1, indices of original Geo polygons that correspond to
%                   each face
%           color:  np x 1 array, where np is the number of polygons in the
%                   corresponding Geo structure, default color index that
%                   will plot each polygon such that it doesn't share a
%                   color with any of its neighbors.

% Copyright 2026 Kelly Kearney

arguments
    opt.reloaddata (1,1) {mustBeNumericOrLogical} =false
    opt.seed (1,1) {mustBeNumeric, mustBeInteger} =123
end

% File where data will be saved

pth = fileparts(mfilename('fullpath'));
suppfol = fullfile(fileparts(pth), "supporting_data");
akdatafile = fullfile(suppfol, 'ak_analysis_polygons.mat');

if ~opt.reloaddata && ~exist(akdatafile, 'file')
    error('Polygon data file does not exist; please run with reloaddata=true');
end

% Calculate or load data

if opt.reloaddata

    % Initialize random number generator to ensure reproducibility

    rng(opt.seed)

    %  EBS/NBS are sampled uniformly by station
    %  GOA/AI/BSS need to be grouped by stratum
    
    Ak = akshapes;
    svys = Ak.base_layers_survey_strata;
    svyg = Ak.base_layers_survey_grid;
    
    % Bering (EBS = 98, NBS = 143) use stations
    % Changes between the 2010 and 2024 EBS grid are very minor; 2024 doesn't
    % include the corner circles so I'm sticking with 2010 here.
    
    Geo.ebs = svyg(svyg.SURVEY_DEFINITION_ID ==  98 & svyg.DESIGN_YEAR == 2010,:);
    Geo.nbs = svyg(svyg.SURVEY_DEFINITION_ID == 143 & svyg.DESIGN_YEAR == 2010,:);
    
    % AI/NBS/GOA use stratified random
    % Note: major overhaul to GOA strata in 2025, using both 2024 and 2025 for
    % now (there is no overlap in stratum number)
    
    Geo.ai =  svys(svys.SURVEY_DEFINITION_ID == 52 & svys.DESIGN_YEAR == 1991,:);
    Geo.goa = svys(svys.SURVEY_DEFINITION_ID == 47 & svys.DESIGN_YEAR == 2024,:);
    Geo.bss = svys(svys.SURVEY_DEFINITION_ID == 78 & svys.DESIGN_YEAR == 2023,:);
    
    % Grids for all regions
    
    Sgrd.ebs = Geo.ebs;
    Sgrd.nbs = Geo.nbs;
    Sgrd.ai =  svyg(svyg.SURVEY_DEFINITION_ID == 52 & svyg.DESIGN_YEAR == 1991,:);
    Sgrd.goa = svyg(svyg.SURVEY_DEFINITION_ID == 47 & svyg.DESIGN_YEAR == 2024,:);
    Sgrd.bss = svyg(svyg.SURVEY_DEFINITION_ID == 78 & svyg.DESIGN_YEAR == 2023,:);
    
    Geo  = structfun(@(x) geotable2struct(x, 'proj', true, 'wrap', 360, 'polyxy', true), Geo,  'uni', 0);
    Sgrd = structfun(@(x) geotable2struct(x, 'proj', true, 'wrap', 360, 'polyxy', true), Sgrd, 'uni', 0);
    
    % Split larger polygons into smaller ones (AI and GOA)
    
    splitreg = ["goa", "ai"];
    
    for ir = 1:length(splitreg)
        rr = splitreg{ir};
    
        Gnew = cell(size(Geo.(rr)));
        nsub = ones(size(Geo.(rr)));
        for ii = 1:length(Geo.(rr))
        
            fprintf('%d/%d\n', ii, length(Geo.(rr)));
        
            sgrdmask = [Sgrd.(rr).AREA_ID] == Geo.(rr)(ii).AREA_ID;
    
            % psub = polysubdivide(Geo.(rr)(ii).polyxy, 2e9, 'nptper', 1000);
            [psub, sgrdidx] = polysubdivide(Geo.(rr)(ii).polyxy, 2e9, 'sgrd', [Sgrd.(rr)(sgrdmask).polyxy]);
    
            nsub(ii) = length(psub);
        
            if nsub(ii) > 1
        
                Gnew{ii} = repmat(Geo.(rr)(ii), length(psub),1);
                for is = 1:nsub(ii)
                    Gnew{ii}(is).Shape = mappolyshape(psub(is).Vertices(:,1), psub(is).Vertices(:,2));
                    Gnew{ii}(is).Shape.ProjectedCRS = Geo.(rr)(ii).Shape.ProjectedCRS;
                    Gnew{ii}(is).AREA_ID = Gnew{ii}(is).AREA_ID + is/100;
                    Gnew{ii}(is).AREA_M2 = area(psub(is));
                    Gnew{ii}(is).Y = [psub(is).Vertices(:,2); NaN]';
                    Gnew{ii}(is).X = [psub(is).Vertices(:,1); NaN]';
                    [Gnew{ii}(is).Lat, Gnew{ii}(is).Lon] = projinv(Geo.(rr)(ii).Shape.ProjectedCRS, Gnew{ii}(is).X, Gnew{ii}(is).Y);
                    Gnew{ii}(is).polyxy = psub(is);
                end

                Tmp = Sgrd.(rr)(sgrdmask);
                for itmp = 1:length(Tmp)
                    Tmp(itmp).AREA_ID = [Geo.(rr)(ii).AREA_ID] + sgrdidx(itmp)/100;
                end
                Sgrd.(rr)(sgrdmask) = Tmp;
    
            end
        end
        Geo.(rr) = [Geo.(rr)(nsub == 1); cat(1, Gnew{:})];
    end
    
    % Patch-plotting triangulation
    
    Tri = structfun(@geostruct2fv, Geo, 'uni', 0);

    % Save

    save(akdatafile, 'Geo', 'Sgrd', 'Tri');

else
    Tmp = load(akdatafile);
    Geo = Tmp.Geo;
    Sgrd = Tmp.Sgrd;
    Tri = Tmp.Tri;
end

end


function T = geostruct2fv(A)

    % Triangulate each polygon

    [f,v] = arrayfun(@(x) poly2fv(x.Lon, x.Lat), A, 'uni', 0);

    % Combine into one f/v patch

    [vunq, ~, iunq] = unique(cat(1, v{:}), 'rows');
    nvertper = cellfun(@(x) size(x,1), v);
    nfaceper = cellfun(@(x) size(x,1), f);
    
    iunq = mat2cell(iunq', 1, nvertper);

    for ii = 1:length(f)
        f{ii} = iunq{ii}(f{ii});
    end
    
    % iscol = cellfun(@(x) isequal(size(x), [3 1]), f); 
    % f(iscol) = cellfun(@(x) x', f(iscol), 'uni', 0);

    % original polygon index for each face

    idx = arrayfun(@(a,b) b*ones(a,1), nfaceper, (1:length(f))', 'uni', 0);

    % Combine into one set of face/vertex/index

    T.faces = cat(1, f{:});
    T.verts = vunq;
    T.fvind = cat(1, idx{:});

    % Default color index: ensures neighboring polygons aren't the same
    % color

    buf = 10; 

    p = [A.polyxy];
    pfat = polybuffer(p, buf);
    np = length(p);

    adj = zeros(np);
    for k = 1:np
        for j = (k+1):np % Only need j for ks we haven't checked yet
            adj(j,k) = area(intersect(pfat(j),pfat(k))); % > 1e-6;
        end
    end

    % [xp,yp] = polycenter(p);
    % [ltp, lnp] = projinv(Geo.(rr)(1).Shape.ProjectedCRS, xp, yp);

    G = graph(adj, 'lower');
    T.color = assigncolors(G, 6);

end

function color = greedycolor(G) 

    n = numnodes(G);
    color = zeros(n,1); % Initialize array to hold color numbering
    numcolors = 1; % Number of colors needed (so far)
    for k = 1:n
        idx = neighbors(G,k); % Get neighbors of the kth node
        idx = idx(idx < k);   % But just those that have an assigned color
        neighborcolors = unique(color(idx)); % Get colors used by neighbors
        % Assign the smallest color value not used by the neighboring nodes
        thiscolor = min(setdiff(1:numcolors,neighborcolors));
        % If there isn't one, add another color to the map
        if isempty(thiscolor)
            numcolors = numcolors + 1;
            thiscolor = numcolors;
        end
        color(k) = thiscolor;
    end
    % disp([num2str(numcolors),' colors needed'])

end

function c = assigncolors(G, maxcol)

    numcolors = Inf;
    n = numnodes(G);
    adj = adjacency(G);
    
    while numcolors > maxcol
        idx = randperm(n);
        Gperm = graph(adj(idx,idx));
        color = greedycolor(Gperm);
        numcolors = max(color);
    end

    c = nan(size(color));
    c(idx) = color;

end

function [psub, sgrdidx] = polysubdivide(pall, maxarea, opt)
% pall: polyshape of larger polygon, can be multifaceted
% maxarea: target area of sub-regions, in square units of p
% nptper: ~number of random points per polygon unit along each axis
%
% This function splits a polyshape object into roughly equally-sized
% regions, following the method outlined here:
% https://blog.cleverelephant.ca/2018/06/polygon-splitting.html 

% TODO: add modified option that builds polygons from subpolys rather than
% voronoi (and seeds points using subpoly centers)

    arguments
        pall
        maxarea (1,1) {mustBeNumeric}
        opt.nptper (1,1) {mustBeNumeric, mustBeInteger} =1000
        opt.sgrd {mustBeA(opt.sgrd, 'polyshape')} =polyshape
    end

    usesgrd = sum([opt.sgrd.NumRegions]) > 0;
    if usesgrd
        opt.sgrd = opt.sgrd(:);
    end

    preg = regions(pall);
    psub = cell(size(preg));

    if usesgrd
         [xgrd, ygrd] = polycenter(opt.sgrd);
    end

    for jj = 1:length(preg)

        % Figure out how many sub-regions we want, based on the target max
        % area

        parea = area(preg(jj));
        npoly = ceil(parea/maxarea);

        if npoly > 1
    
            if usesgrd

                isin = preg(jj).isinterior(xgrd, ygrd);
                xypt = [xgrd(isin,:) ygrd(isin)];

            else
                    
                % Scatter random points within the polygon
    
                xylim = [min(preg(jj).Vertices,[],1); max(preg(jj).Vertices,[],1)];
                dxy = diff(xylim, 1);
                npt = floor(prod(dxy/opt.nptper)); 
                xypt = rand(npt,2).*dxy(1)+xylim(1,:);
                isin = preg(jj).isinterior(xypt);
                xypt = xypt(isin,:);
            end

            % Cluster points into npoly groups

            idx = kmeans(xypt, npoly);

            if usesgrd
                psub{jj} = splitapply(@(x) union(x), opt.sgrd(isin), idx);
            else

                % Find centroid of each group
    
                xyc = splitapply(@(x) mean(x,1), xypt, idx);
    
                % Add extra points far from polygon to allow easier closing of
                % the voronoi regions
    
                prfat = polybuffer(preg(jj), 10e6);
    
                % Calculate voronoi regions around cluster centroids
    
                [v,c] = voronoin([xyc; prfat.Vertices]);
                v(isinf(v)) = 1e9;
                
                % Find intersection of voronoi polygons and the original
                % polygon
                
                psub{jj} = polyshape;
                for ic = size(xyc,1):-1:1
                    psub{jj}(ic) = intersect(polyshape(v(c{ic},1), v(c{ic},2)), preg(jj));
                end
            end

        else
            psub{jj} = preg(jj);
        end

    end

    % Combine all subregions into one polyshape array

    try
        psub = cat(1, psub{:})';
    catch
        psub = cat(2, psub{:})';
    end

    sgrdidx = zeros(size(opt.sgrd));
    for is = 1:length(psub)
        mask = psub(is).isinterior(xgrd, ygrd);
        sgrdidx(mask) = is;
    end

end
