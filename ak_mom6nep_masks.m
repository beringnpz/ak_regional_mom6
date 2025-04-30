% ak_mom6nep_masks

% MOM6-NEP Spring PEEC 2025

%% MOM6-NEP mask variables
%  Build masking variables for MOM6-NEP Alaska-region grid combined with
%  common AFSC polygons.  Eventually, I'd like to move this analysis 
%  into R to streamline use of the akgfmaps package and eliminate the
%  intermediate files, but this will suffice in the meantime.  

%% ... Regions of interest (shapefiles exported from akgfmaps R package)

shapefol = "./akgfmaps_shapefiles";

% ESR regions

esrshape = fullfile(shapefol, 'esr_subareas');
[S.Esr, pcrs] = shapereadprj(esrshape);
S.Esr = wraplon(S.Esr);
for ii = 1:length(S.Esr)
    S.Esr(ii).esrnum = ii;
end

% Bering Sea survey strata

stratashape = fullfile(shapefol, 'surveystrata_ebs');
S.Strata = shapereadprj(stratashape);

% Major survey areas

surveyname = ["ai", "sebs", "nebs", "goa", "ecs"];

S.Survey = arrayfun(@shapereadprj, fullfile(shapefol, "surveyarea_"+surveyname), 'uni', 0);
for ii = 1:length(S.Survey)
    S.Survey{ii}.regionnum = ii;
end

S.Survey = catstruct(1, S.Survey{:});
S.Survey = wraplon(S.Survey);

% Bounds

latlim = minmax([S.Esr.Lat]);
lonlim = minmax(wrapTo360([S.Esr.Lon]));

%% ... quick plot

worldmap(latlim, lonlim)

geoshow(S.Esr, 'facecolor', 'none', 'edgecolor', 'b');
geoshow(S.Survey(1:end-1), 'facecolor', 'none', 'edgecolor', 'r');
geoshow(S.Strata, 'facecolor', 'none', 'edgecolor', 'g', 'linestyle', ':');

%% ... Create Alaska-only cropped version of ocean_static
%
% Building a masking file so I can run some of the basic regional-averaging
% via NCO tools on the gfdl computers.
%
% Static file extracted from Liz's hindcast:
% public5:/home/kak/2025SpringPEEC> tar -xvf /archive/e1n/fre/cefi/NEP/2024_11/NEP_nudge_spinup/gfdl.ncrc6-intel23-repro/history/19930101.nc.tar ./19930101.ocean_daily.nc
%
% then scp'd to klone.  I then extracted the key variables 
%
% ncks -v areacello,deptho,geolat,geolon,wet 19930101.ocean_static.nc ocean_static_alaska.nc

basepath = fullfile(mounteddir('klone'), 'GR011846_reem/kearney/peec2025');
staticfile = fullfile(basepath, '19930101.ocean_static.nc');

Full = buildmasks(staticfile, [], S);

% Start-count-stride to get just Alaska (minus Arctic)
% (Note: I slightly updated the masking to use R akgfmaps package
% boundaries rather than my own collection, which ever so slightly altered
% the sub-region hyperslab.  Setting manually to stay consistent with
% previously-extracted MOM6-NEP hindcast data.
%
% -d xh,1,265 -d yh,442,743 

Scs = struct('xh', [1 265 1], 'yh', [442 302 1]);

% regmask = Full.mask_esrreg <= 3;
% 
% xh = any(regmask,2);
% yh = any(regmask,1);
% Scs = struct('xh', [find(xh,1) sum(xh) 1], ...
%              'yh', [find(yh,1) sum(yh) 1]);

% Repeat with cropped region

Crop = buildmasks(staticfile, Scs, S);

% Command to extract only desired variables from static file

hyperslab = structfun(@(x) sprintf("%d,%d", x(1), x(1)+x(2)-1), Scs, 'uni', 0);
hyperslab =  [fieldnames(hyperslab) struct2cell(hyperslab)]';
hyperslab = sprintf("-d %s,%s ", hyperslab{:});

akstaticfile = fullfile(basepath, "ocean_static_ak.nc");

vstatic = ["areacello", "deptho", "geolat", "geolon", "wet"];
cmd = join(["ncks -F", ...
            "-v " + join(vstatic,","), ...
            hyperslab, ...
            staticfile, ...
            akstaticfile], " ");
system(cmd);

%% ... Add mask variables to file

ncbuild(akstaticfile, Crop.mask_esrreg, ...
    'name', 'mask_esrreg', ...
    'varatts', {'long_name', 'Alaska Ecosystem Status Report (ESR) region', ...
                'flag_values', '1,2,3,4', ...
                'flag_meanings', 'Aleutian Islands, Gulf of Alaska, Bering Sea, Arctic'});

ncbuild(akstaticfile, Crop.mask_survey, ...
    'name', 'mask_survey', ...
    'varatts', {'long_name', 'Alaska groundfish bottom trawl survey region', ...
                'flag_values', '1,2,3,4,5', ...
                'flag_meanings', 'Aleutian Islands (AI), Southeast Bering Sea (SEBS), Northern Bering Sea (NBS), Gulf of Alaska (GOA), Eastern Chukchi Sea (ECS)'});

ncbuild(akstaticfile, Crop.mask_strata, ...
    'name', 'mask_strata', ...
    'varatts', {'long_name', 'Bering Sea groundfish survey strata'});


%% Subs

function [X, p] = shapereadprj(filebase)
    if isstring(filebase); filebase = filebase{1}; end
    p = projcrs(fileread([filebase '.prj']));
    X = shaperead([filebase '.shp']);

    for ii = length(X):-1:1
        [X(ii).Lat, X(ii).Lon] = projinv(p, X(ii).X, X(ii).Y);
    end
end

function S = wraplon(S)
    for ii = 1:length(S)
        S(ii).Lon = wrapTo360(S(ii).Lon);
    end
end

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

    % X.regesr = nan(size(X.Static.geolon));
    % for ii = 1:length(pEsrProj)
    %     tmp = inpolygon(X.Static.geolon, X.Static.geolat, pEsrProj(ii).Vertices(:,2), pEsrProj(ii).Vertices(:,1));
    %     X.regesr(tmp & X.Static.wet==1) = ii;
    % end
    % 
    % % Strata number masks
    % 
    % X.regstrata = nan(size(X.Static.geolon));
    % for ii = 1:length(pStrataProj)
    %     tmp = inpolygon(X.Static.geolon, X.Static.geolat, pStrataProj(ii).Vertices(:,2), pStrataProj(ii).Vertices(:,1));
    %     X.regstrata(tmp & X.Static.wet==1) = stratanum(ii);
    % end
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





