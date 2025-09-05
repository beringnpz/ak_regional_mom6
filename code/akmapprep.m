function [A, Grd] = akmapprep(simname, opt)

arguments
    simname {mustBeTextScalar}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
    opt.maskvar ="mask_esr_area"
    opt.maskval =3
end

% Read relevant static file variables

staticfile = fullfile(opt.datafol, simname, "Level1-2", opt.staticname);

vars = ["geolat", "geolon", "geolat_c", "geolon_c", ...
        "mask_esr_area", "mask_survey_strata", "mask_survey_area", 'wet'];
vars = unique([vars opt.maskvar]);

Grd = ncstruct(staticfile, vars{:});

gsz = size(Grd.geolat);

% EBS region

mask = ismember(Grd.(opt.maskvar), opt.maskval);

% isebs = Grd.mask_esr_area == 3;
latlim = minmax(Grd.geolat(mask), 'expand');
lonlim = minmax(Grd.geolon(mask), 'expand');

% Map decor data: borders, survey region

issebs = Grd.mask_survey_area == 98 & Grd.mask_survey_strata <= 62;
[slon,slat] = mask2poly(Grd.geolon_c, Grd.geolat_c, issebs);

[blat, blon] = deal(cell(1,2));
[blat{1}, blon{1}] = borders('alaska');
[blat{2}, blon{2}] = borders('russia');

warnstate = warning('off', 'MATLAB:polyshape:repairedBySimplify');
bp = polyshape(wrapTo360([blon{:}]), [blat{:}]);
warning(warnstate);

% Projection setup

hfig = figure('visible', 'off');
hb = boxworldmap(latlim, lonlim);
m = getm(hb.ax);
close(hfig);

[xc,yc] = projfwd(m, Grd.geolat_c, Grd.geolon_c);
[bx,by] = projfwd(m, bp.Vertices(:,2), bp.Vertices(:,1));
b = polyshape(bx,by);
[sx, sy] = projfwd(m, slat, slon);

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
A.addboxmap = @() boxworldmap(latlim, lonlim);