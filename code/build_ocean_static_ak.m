function Scs = build_ocean_static_ak(opt)
% BUILD_OCEAN_STATIC_AK Create Alaska-specific version of MOM6 static file
%
% Scs = build_ocean_static_ak(flag);
%
% This script extracts the Alaska subregion of the MOM6-NEP domain from a
% simulation's static file.  It also adds custom variables related to
% Alaska management regions.  See below for the list of new variables and
% their attributes.
%
% Input variables (optional, passed as parameter/value pairs):
%
%   setuponly:  logical scalar, true to run the masking and subsetting
%               calculations but not create the new file
%               [false]
%
%   simname:    name of simulation.  This script expects to find a
%               subfolder with this name under the ak_regional_mom6 data
%               folder.  It will look for the existing static file and
%               create the new Alaska-specific static file in the
%               <datafol>/<simname>/Level1-2/ folder
%               ["mom6nep_hc202507"]
%
%   origname:   base name for original static file, to be appended to
%               <datafol>/<simname>/Level1-2/<simname>_
%               ["ocean_static.nc"]
%
%   newname:    base name for new Alaska-specific static file, to be
%               appended to <datafol>/<simname>/Level1-2/<simname>_
%               ["ocean_static_ak.nc"]
%
%   datafol:    CEFI data folder path.  Default is the path returned by the
%               cefidatafol.m function
%
% Output variables:
%
%   Scs:        ncstruct subsetting structure, holding [start count stride]
%               indices used to subset each dimension in the original
%               static file for the Alaska region.
%
% New variables added to copied-and-trimmed version of the static file:
%
%   mask_adfg_stat_area    
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'ADFG Stat Area'
%   mask_bsierp_region     
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Marine Regions of the Bering Sea Shelf and Slope (Bering Sea Integrated Ecosystem Research Program), https://doi.org/10.5065/D6DF6P6C'
%                      flag_values   = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16'
%                      flag_meanings = '1) AK peninsula, 
%                                       2) South inner shelf, 
%                                       3) South middle shelf, 
%                                       4) South outer shelf, 
%                                       5) Pribilofs, 
%                                       6) Central middle shelf, 
%                                       7) Central inner shelf, 
%                                       8) North outer shelf, 
%                                       9) St. Matthew, 
%                                       10) North middle shelf, 
%                                       11) North inner shelf, 
%                                       12) St. Lawrence, 
%                                       13) South Bering Strait, 
%                                       14) Norton Sound,  
%                                       15) Off-shelf north, 
%                                       16) Off-shelf southeast'
%   mask_crab_strata_BKC   
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Blue King Crab survey strata'
%                      flag_values   = '1,2,3,4'
%                      flag_meanings = '1) Pribilof MTCA (Pribilof Islands BKC), 
%                                       2) Pribilof Single (Pribilof Islands BKC), 
%                                       3) St. Matthew Single (St. Matthew BKC), 
%                                       4) St. Matthew MTCA (St. Matthew BKC)'
%   mask_crab_strata_RKC   
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Red King Crab survey strata'
%                      flag_values   = '1,2,3,4'
%                      flag_meanings = '1) Bristol Bay (Bristol Bay RKC), 
%                                       2) Pribilof MTCA (Pribilof Islands RKC), 
%                                       3) Pribilof Single (Pribilof Islands RKC), 
%                                       4) Norton Sound (Norton Sound RKC)'
%   mask_crab_strata_snow  
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Snow crab survey strata'
%                      flag_values   = '1,2,3'
%                      flag_meanings = '1) Single (EBS snow crab), 
%                                       2) Pribilof MTCA (EBS snow crab), 
%                                       3) St. Matthew MTCA (EBS snow crab)'
%   mask_crab_strata_tanner
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Tanner crab survey strata'
%                      flag_values   = '1,2,3,4'
%                      flag_meanings = '1) East 166 (Tanner E), 
%                                       2) West 166 (Tanner W), 
%                                       3) Pribilof MTCA (Tanner W), 
%                                       4) St. Matthew MTCA (Tanner W)'
%   mask_esr_area          
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Alaska Marine Ecosystem Status Reports: Ecosystem Area'
%                      flag_values   = '1,2,3,4'
%                      flag_meanings = '1) Aleutian Islands, 
%                                       2) Arctic, 
%                                       3) Eastern Bering Sea, 
%                                       4) Gulf of Alaska'
%   mask_esr_subarea       
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Alaska Marine Ecosystem Status Reports: Ecosystem Subarea'
%                      flag_values   = '1,2,3,4,5,6,7,8,9'
%                      flag_meanings = '1) Central Aleutians, 
%                                       2) Western Aleutians, 
%                                       3) Northern Bering Sea, 
%                                       4) Southeastern Bering Sea, 
%                                       5) Eastern Gulf of Alaska, 
%                                       6) Western Gulf of Alaska, 
%                                       7) Eastern Aleutians, 
%                                       8) Western Gulf of Alaska (inside), 
%                                       9) Eastern Gulf of Alaska (inside)'
%   mask_inpfc_strata      
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'INPFC Area'
%                      flag_values   = '929,939,919,959,949,6499,5699,799,299'
%                      flag_meanings = '929) Chirikof, 
%                                       939) Kodiak, 
%                                       919) Shumagin, 
%                                       959) Southeastern, 
%                                       949) Yakutat, 
%                                       6499) Central Aleutians, 
%                                       5699) Eastern Aleutians, 
%                                       799) Southern Bering Sea, 
%                                       299) Western Aleutians'
%   
%   mask_nmfs_area         
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'NMFS Refixporting Area (Electronic Code of Federal Regulations, Part 679, Appendix A)'
%                      flag_values   = '400,508,509,512,513,514,516,517,518,519,521,523,524,530,541,542,543,550,610,620,630,640,649,650,659'
%                      flag_meanings = '400) Chukchi Sea, 
%                                       518) Bogoslof District,
%                                       541) Eastern Aleutian District, 
%                                       542) Central Aleutian District, 
%                                       543) Western Aleutian District, 
%                                       550) Donut Hole, 
%                                       610) Western GOA, Shumagin District, 
%                                       620) Central GOA, Chirkof District, 
%                                       630) Central GOA, Kodiak District, 
%                                       640) Eastern GOA, West Yakutat District, 
%                                       649) Prince William Sound, 
%                                       650) Eastern GOA, Southeast Outside District, 
%                                       659) Eastern GOA, Southeast Inside District'
%   mask_survey_area       
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Bottom trawl survey area'
%                      flag_values   = '6,47,52,78,98,143'
%                      flag_meanings = '6) Chukchi Sea Trawl Survey (2012), 
%                                       47) Gulf of Alaska Bottom Trawl Survey (2025), 
%                                       52) Aleutian Islands Bottom Trawl Survey (1991), 
%                                       78) Eastern Bering Sea Slope Bottom Trawl Survey (2023), 
%                                       98) Eastern Bering Crab/Groundfish Bottom Trawl Survey (2022), 
%                                       143) Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension (2022)'
%   mask_survey_strata     
%          Size:       342x297
%          Dimensions: ih,jh
%          Datatype:   single
%          Attributes:
%                      long_name     = 'Bottom trawl survey strata ID'

% Copyright 2025 Kelly Kearney

arguments
    opt.setuponly (1,1) {mustBeNumericOrLogical} =false
    opt.simname  {mustBeTextScalar} ="mom6nep_hc202507"
    opt.origname {mustBeTextScalar} ="ocean_static.nc"
    opt.newname  {mustBeTextScalar} ="ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
end

if nargin < 1
    flag = false;
end
validateattributes(flag, {'logical'}, {'scalar'});

%--------------------
% Folder setup
%--------------------

% Data location

origfile = fullfile(opt.datafol, opt.simname, "Level1-2", opt.simname+"_"+opt.origname);
newfile  = fullfile(opt.datafol, opt.simname, "Level1-2", opt.simname+"_"+opt.newname);

%--------------------
% Read data
%--------------------

% Regions of interest (from akgfmaps R package)

S = akshapes;

% ... Bounds

latlim = minmax(S.esr_regions.Shape.InternalData.VertexCoordinate1, 'expand');
lonlim = minmax(wrapTo360(S.esr_regions.Shape.InternalData.VertexCoordinate2), 'expand');

% Static file coordinates

Grd = ncstruct(origfile, 'geolat', 'geolon', 'geolat_c', 'geolon_c');

usept = ~isnan(Grd.geolat) & ...
        Grd.geolat >= latlim(1) & Grd.geolat <= latlim(2) & ...
        Grd.geolon >= lonlim(1) & Grd.geolon <= lonlim(2);

staticg = geopointshape(Grd.geolat(usept), Grd.geolon(usept));

pcrs = projcrs(3338);
[x,y] = projfwd(pcrs, Grd.geolat(usept), Grd.geolon(usept));
staticm = mappointshape(x, y);
staticm.ProjectedCRS = pcrs;

%--------------------
% Masking
%--------------------

sz = size(Grd.geolat);

% ESR regions

Mask.esr_area.val = val2mask(interpgeotable(S.esr_regions(1:4,:), staticg, [], 0.001), sz, usept, true);
Mask.esr_area.long_name = "Alaska Marine Ecosystem Status Reports: Ecosystem Area";
Mask.esr_area.flag_meanings = S.esr_regions.Area_Name(1:4);

% ESR subregions

Mask.esr_subarea.val = val2mask(interpgeotable(S.esr_regions(5:11,:), staticg, [], 0.001), sz, usept, true);
tmp = val2mask(interpgeotable(S.esr_regions(13:14,:), staticg, [], 0.001), sz, usept, true);

isin = ~isnan(tmp) & isnan(Mask.esr_subarea.val);
Mask.esr_subarea.val(isin) = tmp(isin) + 7;

Mask.esr_subarea.long_name = "Alaska Marine Ecosystem Status Reports: Ecosystem Subarea";
Mask.esr_subarea.flag_meanings = [S.esr_regions.Area_Name(5:11); S.esr_regions.Area_Name(13:14)+" (inside)"];

% ADFG areas

Mask.adfg_stat_area.val = val2mask(interpgeotable(S.adfg_areas, staticg, 'AreaID'), sz, usept, true);
Mask.adfg_stat_area.long_name = "ADFG Stat Area";

% Survey area (most recent design year)

S.base_layers_survey_area = sortrows(S.base_layers_survey_area, {'SURVEY_DEF', 'DESIGN_YEA'});
[~,iunq] = unique(S.base_layers_survey_area.SURVEY_DEF, 'last');

Mask.survey_area.val = val2mask(interpgeotable(S.base_layers_survey_area(iunq,:), staticm, 'SURVEY_DEF', 100), sz, usept, true);
Mask.survey_area.long_name = "Bottom trawl survey area";
Mask.survey_area.flag_values = S.base_layers_survey_area.SURVEY_DEF(iunq);
Mask.survey_area.flag_meanings = S.base_layers_survey_area.SURVEY_NAM(iunq) + compose(" (%d)", S.base_layers_survey_area.DESIGN_YEA(iunq));

% Survey strata (most recent design year)

isin = ismember(S.base_layers_survey_strata{:,{'SURVEY_DEF','DESIGN_YEA'}}, S.base_layers_survey_area{iunq,{'SURVEY_DEF','DESIGN_YEA'}}, 'rows');

Mask.survey_strata.val = val2mask(interpgeotable(S.base_layers_survey_strata(isin,:), staticm, "AREA_ID", 100), sz, usept, true);
Mask.survey_strata.long_name = "Bottom trawl survey strata ID";

% INPFC strata 

Mask.inpfc_strata.val = val2mask(interpgeotable(S.base_layers_inpfc_strata, staticm, 'AREA_ID'), sz, usept, true);
Mask.inpfc_strata.long_name = "INPFC Area";
Mask.inpfc_strata.flag_values = S.base_layers_inpfc_strata.AREA_ID;
Mask.inpfc_strata.flag_meanings = S.base_layers_inpfc_strata.AREA_NAME;

% BSIERP regions

S.bsierp_regions = sortrows(S.bsierp_regions, 'AreaID');

Mask.bsierp_region.val = val2mask(interpgeotable(S.bsierp_regions, staticg, 'AreaID', 0.001), sz, usept, true);
Mask.bsierp_region.long_name = "Marine Regions of the Bering Sea Shelf and Slope (Bering Sea Integrated Ecosystem Research Program), https://doi.org/10.5065/D6DF6P6C";
Mask.bsierp_region.flag_values = S.bsierp_regions.AreaID;
Mask.bsierp_region.flag_meanings = S.bsierp_regions.Area_Name;

% Crab strata

Mask.crab_strata_RKC.val = val2mask(interpgeotable(S.crab_strata_RKC, staticm, [], 100), sz, usept, true);
Mask.crab_strata_RKC.long_name = "Red King Crab survey strata";
Mask.crab_strata_RKC.flag_meanings = S.crab_strata_RKC.STRATUM + " (" +  S.crab_strata_RKC.STOCK + ")";

Mask.crab_strata_BKC.val = val2mask(interpgeotable(S.crab_strata_BKC, staticm, [], 100), sz, usept, true);
Mask.crab_strata_BKC.long_name = "Blue King Crab survey strata";
Mask.crab_strata_BKC.flag_meanings = S.crab_strata_BKC.STRATUM + " (" +  S.crab_strata_BKC.STOCK + ")";

Mask.crab_strata_snow.val = val2mask(interpgeotable(S.crab_strata_snow, staticm, [], 100), sz, usept, true);
Mask.crab_strata_snow.long_name = "Snow crab survey strata";
Mask.crab_strata_snow.flag_meanings = S.crab_strata_snow.STRATUM + " (" +  S.crab_strata_snow.STOCK + ")";

Mask.crab_strata_tanner.val = val2mask(interpgeotable(S.crab_strata_tanner, staticm, [], 100), sz, usept, true);
Mask.crab_strata_tanner.long_name = "Tanner crab survey strata";
Mask.crab_strata_tanner.flag_meanings = S.crab_strata_tanner.STRATUM + " (" +  S.crab_strata_tanner.STOCK + ")";

% NMFS Areas

Mask.nmfs_area.val = val2mask(interpgeotable(S.nmfs_areas, staticg, 'AreaID', 0.001), sz, usept, true);
Mask.nmfs_area.long_name = "NMFS Reporting Area (Electronic Code of Federal Regulations, Part 679, Appendix A)";
Mask.nmfs_area.flag_values = S.nmfs_areas.AreaID;
Mask.nmfs_area.flag_meanings = S.nmfs_areas.Area_Name;

%--------------------
% Test plot
%--------------------

fld = fieldnames(Mask);
for ii = 1:length(fld)
    figure('name', fld{ii});
    worldmap(latlim, lonlim)
    [ic, cval] = findgroups(Mask.(fld{ii}).val(:));
    ic = reshape(ic, size(Mask.(fld{ii}).val));
    nc = max(ic(:));

    pcolorm(Grd.geolat_c, Grd.geolon_c, padend(ic));
    cb = colorbar;
    set(cb, 'ticks', 1:nc, 'ticklabels', compose("%g", cval));
    set(gca, 'clim', [0.5 nc+0.5], 'colormap', distinguishable_colors(nc));
end

%--------------------
% Create extended and
% cropped static file
%--------------------

% Add mask variables

if ~opt.setuponly

    copyfile(origfile, newfile);

    fld = fieldnames(Mask);
    for ii = 1:length(fld)

        if isfield(Mask.(fld{ii}), 'flag_meanings')

            nflag = length(Mask.(fld{ii}).flag_meanings);
            if isfield(Mask.(fld{ii}), 'flag_values')
                flagval = Mask.(fld{ii}).flag_values;
            else
                flagval = (1:nflag)';
            end
               
            flagstr = compose("%d) ", flagval) + Mask.(fld{ii}).flag_meanings;

            ncbuild(newfile, single(Mask.(fld{ii}).val), ...
                'name', "mask_" + fld{ii}, ...
                'varatts', {'long_name', Mask.(fld{ii}).long_name, ...
                            'flag_values',   join(compose("%d", flagval),","), ...
                            'flag_meanings', join(flagstr,                ", ")});
        else

            ncbuild(newfile, single(Mask.(fld{ii}).val), ...
                'name', "mask_" + fld{ii}, ...
                'varatts', {'long_name', Mask.(fld{ii}).long_name});

        end
    end

end

% Extract relevant hyperslab only

xh = any(~isnan(Mask.esr_area.val),2);
yh = any(~isnan(Mask.esr_area.val),1);

% start-count-stride

Scs = struct('ih', [find(xh,1) sum(xh) 1], ...
             'jh', [find(yh,1) sum(yh) 1]);
Scs.iq = Scs.ih + [0 1 0];
Scs.jq = Scs.jh + [0 1 0];

% Extract subset with NCKS

cmd = join(["ncks -F -O", ...
            hyperslabstring(Scs), ...
            newfile, ...
            newfile], " ");

if ~opt.setuponly
    system(cmd);
end


end


% ******** Subfunctions ********

function varargout = interpgeotable(gtbl, pts, prop, dx)
%INTERPGEOTABLE Interpolate point values based on which polygons are in
%
% Input variables:
%
%   gtbl:   geospatial table where Shape column holds geopolyshapes or
%           mappolyshapes
%
%   pts:    geopointshape or mappointshape, with conforming projection
%           details to the gtbl input
%
%   prop:   string array, table column name(s) to interpolate to each point
%           location.  Values of this column must be numeric scalars. 
%
%   dx:     (optional) tolerance used to simplify polygon vertices via the
%           Douglas-Peuker polyline simplification algorithm.
%           Simplification can greatly reduce in-polygon calculations time
%           when polygon vertex resolution is much higher than necessary
%           for the application.  If empty or not included, no
%           simplification will be performed.
%
% Output variables:
%
%   numerical array of same size as pts, holding value prop value
%   corresponding to each point location.

    % Input checks

    isgeo = isa(gtbl.Shape, 'geopolyshape') && isa(pts, 'geopointshape');
    ismap = isa(gtbl.Shape, 'mappolyshape') && isa(pts, 'mappointshape');

    if ischar(prop)
        prop = string(prop);
    end

    % Convert to geo/mapstruct (allows for line simplification and more
    % efficient inpolygon calculations)

    S = geotable2struct(gtbl);

    % Check if points are in each polygon

    if isgeo
        [lt,ln] = vertexdata(pts);
        ln = wrapTo360(ln);

        if nargin > 3
            for is = 1:length(S)
                ps = dpsimplify([S(is).Lon', S(is).Lat'], dx);
                S(is).Lon = ps(:,1)';
                S(is).Lat = ps(:,2)';
            end
        end

        idx = zeros(size(lt));
        cnt = zeros(size(lt));

        for ii = 1:length(S)
            tmp = inpolygon(ln, lt, S(ii).Lon, S(ii).Lat);
            idx(tmp==1) = ii;
            cnt(tmp==1) = cnt(tmp==1)+1;
        end

    elseif ismap
        [x,y] = vertexdata(pts);

        if nargin > 3
            for is = 1:length(S)
                ps = dpsimplify([S(is).X', S(is).Y'], dx);
                S(is).X = ps(:,1)';
                S(is).Y = ps(:,2)';
            end
        end

        idx = zeros(size(x));
        cnt = zeros(size(y));

        for ii = 1:length(S)
            tmp = inpolygon(x, y, S(ii).X, S(ii).Y);
            idx(tmp==1) = ii;
            cnt(tmp==1) = cnt(tmp==1)+1;
        end
    else
        error('Unexpected data types');
    end

    % Any overlaps?

    if any(cnt(:)>1)
        warning('Overlapping polygons found; assigning values from last polygon checked')
    end

    % Assign properties

    if nargin < 3 || isempty(prop)
        varargout{1} = idx;
    else
        mask = idx > 0;
        for ip = 1:length(prop)
            pval = zeros(size(idx));
            pval(mask) = gtbl.(prop{ip})(idx(mask));
            varargout{ip} = pval;
        end
    end
end

function [x,y] = vertexdata(shp)
%VERTEXDATA Extract coordinates from geoshape or mapshape objects 

    x = shp.InternalData.VertexCoordinate1;
    y = shp.InternalData.VertexCoordinate2;

end

function a = geotable2struct(b)
%VERTEXDATA Convert geospatial table to a geographic data structure

    a = table2struct(b);

    for ii = 1:length(a)
        x = a(ii).Shape.InternalData.VertexCoordinate1;
        y = a(ii).Shape.InternalData.VertexCoordinate2;
        idx = a(ii).Shape.InternalData.IndexOfLastVertex;

        n = idx - [1 idx(1:end-1)+1]+1;

        x = mat2cell(x, 1, n);
        y = mat2cell(y, 1, n);
        x = cellfun(@(z) [z NaN], x, 'uni', 0);
        y = cellfun(@(z) [z NaN], y, 'uni', 0);
        x = cat(2, x{:});
        y = cat(2, y{:});

        if isa(a(ii).Shape, 'geopolyshape')
            a(ii).Lat = x;
            a(ii).Lon = wrapTo360(y);
        elseif isa(a(ii).Shape, 'mappolyshape')
            a(ii).Y = y;
            a(ii).X = x;
        end
    end

end

function y = val2mask(x, sz, usept, replace0)
%VAL2MASK Place vector of values in an array based on a logical mask

    y = zeros(sz);
    y(usept) = x;
    if replace0
        y(y == 0) = NaN;
    end

end

function x = padend(x)
%PADEND Add trailing row and column of NaNs to a 2D array

    x = [x nan(size(x,1),1); nan(1,size(x,2)+1)];

end

function str = hyperslabstring(Scs)
%HYPERSLABSTRING Create NCO-style hyperslab string based on ncstruct
%start-count-stride structure

    Hs = structfun(@(x) sprintf("%d,%d",[x(1) x(1)+x(2)-1]), Scs, 'uni', 0);
    str =  [fieldnames(Hs) struct2cell(Hs)]';
    str = sprintf("-d %s,%s ", str{:});
end





