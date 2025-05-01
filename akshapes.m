function [S, pcrs] = akshapes(varargin)
%AKSHAPES Simple wrapper to read in akgfmaps shapefiles
%
% S = akshapes
%
% This function reads in the shapefiles exported from the akgfmaps R
% library (see export_akgfmaps_polygons.R).  
%
% Output variables:
%
%   S:  structure, fields are each geographic data structures corresponding
%       to various akgfmaps features:
%
%       Esr:    Ecosystem Status Report (ESR) regions
%
%       Strata: Bering Sea groundfish survey strata
%
%       Survey: groundfish survey regions

% Copyright 2025 Kelly Kearney

pth = fileparts(mfilename('fullpath'));
shapefol = fullfile(pth, "akgfmaps_shapefiles");

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

end

%-------------
% Subfunctions
%-------------

% Read projected shapefile

function [X, p] = shapereadprj(filebase)
    if isstring(filebase); filebase = filebase{1}; end
    p = projcrs(fileread([filebase '.prj']));
    X = shaperead([filebase '.shp']);

    for ii = length(X):-1:1
        [X(ii).Lat, X(ii).Lon] = projinv(p, X(ii).X, X(ii).Y);
    end
end

% Wrap logitude to 0-360

function S = wraplon(S)
    for ii = 1:length(S)
        S(ii).Lon = wrapTo360(S(ii).Lon);
    end
end

% Create polygons

function ps = shape2poly(S, pcrs, m)

    for ii = length(S):-1:1
        psxy(ii) = polyshape(S(ii).X, S(ii).Y);
    end    
    

end

