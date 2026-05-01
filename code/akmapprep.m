function [A, Grd] = akmapprep(opt)
%AKMAPPREP Calulate map-related parameters
%
% [A, Grd] = akmapprep(opt)
%
% This function reads in grid information and does some preliminary
% calculations to assist in plotting MOM6-based maps for the Alaska region.
%
% Optional input variables:
%
%   cpopts:     cefi portal options object corresponding to the simulation
%               
%   gfiles:     1 x n string array of filenames where relevent grid and
%               masking data can be found.  If not specified, it will use
%               the grid=static, variable=ocean_static and grid=extra,
%               variable=ak_masks files from the specified simulation.
%
%   maskvar:    masking variable (in static collection) used to set map
%               latitude and longitude limits.  The range will include all
%               grid cell where maskvar == maskval.
%               ["mask_esr_area"]
%
%   maskval:    masking variable value used to set map latitude and
%               longitude limits. The range will include all grid cell
%               where maskvar matches any maskval. 
%               [3]
%
% Output variables:
%
%   A:          1x1 structure holding the following fields:
%
%               xc:     x-coordinates of grid corner points, projected
%                       using the default worldmap option for the indicated
%                       lat/lon limits
%
%               yc:     y-coordinates of grid corner points, projected
%                       using the default worldmap option for the indicated
%                       lat/lon limits
%
%               bx:     x-coordinate of projected coastline border polygon
%
%               by:     y-coordinate of projected coastline border polygon
%
%               b:      projected border polygon polyshape objected
%
%               sx:     x-coordinate of projected EBS survey region polygon
%
%               sy:     y-coordinate of projected EBS survey region polygon
%
%               latlim: 1x2 array, latitude limits of selected mask
%
%               lonlim: 1x2 array, longitude limits of selected mask
%
%               addboxmap: function handle to create a boxed-out worldmap
%                       (see boxworldmap) encompassing the mask region
%
%   Grd:        1x1 structure with grid-related static variables
%
%               geolat:             latitude of middle h-points
%
%               geolat_c:           latitude of corner q-points
%
%               geolon:             longitude of middle h-points
%
%               geolon_c:           longitude of corner q-points
%
%               mask_esr_area:      ESR areas: 
%                                   1) Aleutian Islands, 
%                                   2) Arctic, 
%                                   3) Eastern Bering Sea, 
%                                   4) Gulf of Alaska'
%
%               mask_survey_area:   Survey areas:
%                                   6)   Chukchi Sea Trawl Survey (2012), 
%                                   47)  Gulf of Alaska Bottom Trawl Survey
%                                        (2025) 
%                                   52)  Aleutian Islands Bottom Trawl
%                                        Survey (1991) 
%                                   78)  Eastern Bering Sea Slope Bottom
%                                        Trawl Survey (2023) 
%                                   98)  Eastern Bering Crab/Groundfish
%                                        Bottom Trawl Survey (2022) 
%                                   143) Northern Bering Sea
%                                        Crab/Groundfish Survey - Eastern
%                                        Bering Sea Shelf Survey Extension
%                                        (2022)'
%
%               mask_survey_strata: Bottom trawl survey strata ID
%
%               wet:                ocean mask (1=water, 0=land)

% Copyright 2026 Kelly Kearney

arguments
    opt.cpopts (1,1) {mustBeA(opt.cpopts, "cefiportalopts")} =cefiportalopts()
    opt.gfiles {mustBeText} =''
    opt.maskvar ="mask_esr_area"
    opt.maskval =3
    opt.latbuffer = 0.1;
    opt.lonbuffer = 0.1;
end

% Read relevant static file variables

vars = ["geolat", "geolon", "geolat_c", "geolon_c", ...
        "mask_esr_area", "mask_survey_strata", "mask_survey_area", ...
        "wet"];
vars = unique([vars opt.maskvar]);

Grd = readcefigridvars(opt.cpopts, vars);

% Region of interest

mask = ismember(Grd.(opt.maskvar), opt.maskval);

[latlim, lonlim] = deal(nan(1,2));

[latlim(1),latlim(2)] = bounds(Grd.geolat(mask));
[lonlim(1),lonlim(2)] = bounds(Grd.geolon(mask));

latlim = latlim + [-1 1]*opt.latbuffer*diff(latlim);
lonlim = lonlim + [-1 1]*opt.lonbuffer*diff(lonlim);
% latlim = minmax(Grd.geolat(mask), 'expand');
% lonlim = minmax(Grd.geolon(mask), 'expand');

% Map decor data: borders, survey region

issebs = Grd.mask_survey_area == 98 & Grd.mask_survey_strata <= 62;
[slon,slat] = mask2poly(Grd.geolon_c, Grd.geolat_c, issebs);

[blat, blon] = deal(cell(1,3));
[blat{1}, blon{1}] = borders('alaska');
[blat{2}, blon{2}] = borders('russia');
[blat{3}, blon{3}] = borders('canada');

warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
bp = polyshape(wrapTo360([blon{:}]), [blat{:}]);

% Projection setup

warning('off', 'map:projections:notStandardProjection');
hfig = figure('visible', 'off');
hb = boxworldmap(latlim, lonlim);
m = getm(hb.ax(1));
close(hfig);

[xc,yc] = projfwd(m, Grd.geolat_c, Grd.geolon_c);
[bx,by] = projfwd(m, bp.Vertices(:,2), bp.Vertices(:,1));
b = polyshape(bx,by);
[sx, sy] = projfwd(m, slat, slon);

warning(warnstate);

% Output

A.xc = xc;
A.yc = yc;
A.bx = bx;
A.by = by;
A.sx = sx;
A.sy = sy;
A.b = b;
A.latlim = latlim;
A.lonlim = lonlim;
A.addboxmap = @(x) boxworldmap(latlim, lonlim, x{:});