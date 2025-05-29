function Scs = build_ocean_static_ak(flag)
% BUILD_OCEAN_STATIC_AK Create Alaska-specific version of MOM6 static file
%
% Scs = build_ocean_static_ak(flag);
%
% This script extracts the Alaska subregion of the MOM6-NEP domain from a
% simulation's static file.  It also adds a few custom variables related to
% Alaska management regions:
%
%   mask_esrreg
%          Size:       254x297
%          Dimensions: xh,yh
%          Datatype:   double
%          Attributes:
%                      long_name     = 'Alaska Ecosystem Status Report (ESR) region'
%                      flag_values   = '1,2,3,4,5,6,7'
%                      flag_meanings = 'Central Aleutians, 
%                                       Western Aleutians, 
%                                       Northern Bering Sea, 
%                                       Southeastern Bering Sea, 
%                                       Eastern Gulf of Alaska, 
%                                       Western Gulf of Alaska, 
%                                       Eastern Aleutians'
%   mask_survey
%          Size:       254x297
%          Dimensions: xh,yh
%          Datatype:   double
%          Attributes:
%                      long_name     = 'Alaska groundfish bottom trawl survey region'
%                      flag_values   = '1,2,3,4,5'
%                      flag_meanings = 'Aleutian Islands (AI), 
%                                       Southeast Bering Sea (SEBS), 
%                                       Northern Bering Sea (NBS), 
%                                       Gulf of Alaska (GOA), 
%                                       Eastern Chukchi Sea (ECS)'
%   mask_strata
%          Size:       254x297
%          Dimensions: xh,yh
%          Datatype:   double
%          Attributes:
%                      long_name = 'Bering Sea groundfish survey strata'
%
% Note: I had originally hoped to complete this workflow in R to streamline
% the use of the AFSC GAP program's akgfmaps package (which standardizes
% distribution of these management regions), eliminate intermediate files, 
% and allow for easy updates if the boundaries change or we want to add a
% new one.  The tight timeline (and my minimal R knowledge) prevented this
% for now.
%
% Input variables:
%
%   flag:   true to create new file, false to return the subsetting info
%           without creating the file.  Default is false

if nargin < 1
    flag = false;
end
validateattributes(flag, {'logical'}, {'scalar'});

% Regions of interest (shapefiles exported from akgfmaps R package)

[S, pcrs] = akshapes;

% Bounds

latlim = minmax([S.Esr.Lat]);
lonlim = minmax(wrapTo360([S.Esr.Lon]));

% Data location

datafolfile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'simulation_data', 'data_folder.txt');
if ~exist(datafolfile, 'file')
    error('No data_folder.txt file found');
end
datafol = fileread(datafolfile);

simname = 'mom6nep_hc202411';
staticfile = fullfile(datafol, simname, 'Level0', '19930101.ocean_static.nc');
lev1 = fullfile(datafol, simname, 'Level1-2');

% Create masks

Full = buildmasks(staticfile, [], S);

% Start-count-stride to get just Alaska (minus Arctic)
% Scs = struct('xh', [1 265 1], 'yh', [442 302 1]);

regmask = ~isnan(Full.mask_esrreg);

xh = any(regmask,2);
yh = any(regmask,1);
Scs = struct('xh', [find(xh,1) sum(xh) 1], ...
             'yh', [find(yh,1) sum(yh) 1]);
Scs.xq = Scs.xh + [0 1 0];
Scs.yq = Scs.yh + [0 1 0];

% Repeat with cropped region

Crop = buildmasks(staticfile, Scs, S);

% Command to extract only desired variables from static file

hyperslab = structfun(@(x) sprintf("%d,%d", x(1), x(1)+x(2)-1), Scs, 'uni', 0);
hyperslab =  [fieldnames(hyperslab) struct2cell(hyperslab)]';
hyperslab = sprintf("-d %s,%s ", hyperslab{:});

akstaticfile = fullfile(lev1, "ocean_static_ak.nc");

vstatic = ["areacello", "deptho", "geolat", "geolon", "wet", "geolon_c", "geolat_c"];
cmd = join(["ncks -F", ...
            "-v " + join(vstatic,","), ...
            hyperslab, ...
            staticfile, ...
            akstaticfile], " ");
if flag
    system(cmd);
end

% Add mask variables to file

if flag

    ncbuild(akstaticfile, Crop.mask_esrreg, ...
        'name', 'mask_esrreg', ...
        'varatts', {'long_name', 'Alaska Ecosystem Status Report (ESR) region', ...
                    'flag_values', join(compose("%d", 1:length(S.Esr)),","), ...
                    'flag_meanings', join(string({S.Esr.AREA_NAME}),", ")});
    
    ncbuild(akstaticfile, Crop.mask_survey, ...
        'name', 'mask_survey', ...
        'varatts', {'long_name', 'Alaska groundfish bottom trawl survey region', ...
                    'flag_values', '1,2,3,4,5', ...
                    'flag_meanings', 'Aleutian Islands (AI), Southeast Bering Sea (SEBS), Northern Bering Sea (NBS), Gulf of Alaska (GOA), Eastern Chukchi Sea (ECS)'});
    
    ncbuild(akstaticfile, Crop.mask_strata, ...
        'name', 'mask_strata', ...
        'varatts', {'long_name', 'Bering Sea groundfish survey strata'});

end

end


%% Subs

function X = buildmasks(staticfile, Scs, S)

    if isempty(Scs)
        X.Static = ncstruct(staticfile, 'wet', 'deptho', 'geolat', 'geolon', 'xh', 'yh');
    else
        X.Static = ncstruct(staticfile, 'wet', 'deptho', 'geolat', 'geolon', 'xh', 'yh', Scs);
    end
    
    X.latlim = minmax(X.Static.geolat);
    X.lonlim = minmax(X.Static.geolon);
    
    % ESR region masks (1 = Aleutians, 2 = GoA, 3 = Bering, 4 = Arctic)
    
    X.mask_esrreg = interpshapefile(S.Esr,    X.Static.geolat, X.Static.geolon, 'esrnum');
    X.mask_strata = interpshapefile(S.Strata, X.Static.geolat, X.Static.geolon, 'STRATUM');
    X.mask_survey = interpshapefile(S.Survey, X.Static.geolat, X.Static.geolon, 'regionnum');

end

function pshape = shape2poly(S,p)
    if nargin == 1 && isfile([S '.shp']) && isfile([S '.prj'])
        p = projcrs(fileread([S '.prj']));
        S = shaperead([S '.shp']);
    end

    for ii = length(S):-1:1
        [S(ii).Lat, S(ii).Lon] = projinv(p, S(ii).X, S(ii).Y);
        pshape(ii) = polyshape(S(ii).Lat, wrapTo360(S(ii).Lon));
    end
end





