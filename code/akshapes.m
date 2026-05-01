function S = akshapes(opt)
%AKSHAPES Simple wrapper to read in akgfmaps geometries
%
% S = akshapes
%
% This function reads in the polygon geometries from the akgfmaps R
% library.  It reads data directly where possible, relying on a few
% intermediate shapefiles where necessary.
%
% Input variables (passed as parameter/value pairs):
%
%   akgfmapspath:   path to local cloned copy of akgfmaps repository
%                   (git@github.com:afsc-gap-products/akgfmaps.git)
%                   ["/Users/kelly.kearney/Documents/Repos/Other/akgfmaps/"]
%
%   reloaddata:     logical scalar, true to reread data from file and
%                   recreate the locally-saved .mat file, false to read
%                   .mat file.
%                   [false]
%
%   dataset:        string array, names of datasets to import.  If empty,
%                   all available tables will be imported.  See S
%                   description for options.
%
% Output variables:
%
%   S:          structure of geospatial tables corresponding to the
%               datasets. Some tables include geopolyshapes while others
%               have mappolyshapes; this combination preserves the
%               projection data (or lack thereof) in the original files.
%               Currently, this structure includes the following fields:
% 
%                              adfg_areas: [1736×40 table]
%               base_layers_survey_strata: [331×6 table]
%                 base_layers_survey_area: [11×7 table]
%                base_layers_inpfc_strata: [9×7 table]
%                    base_layers_longline: [311×7 table]
%                 base_layers_ll_stations: [103×17 table]
%                          bsierp_regions: [17×40 table]
%                         crab_strata_RKC: [4×4 table]
%                         crab_strata_BKC: [4×4 table]
%                      crab_strata_tanner: [4×4 table]
%                        crab_strata_snow: [3×4 table]
%                             esr_regions: [14×40 table]
%                        inpfc_goa_strata: [59×3 table]
%                         inpfc_ai_strata: [864×7 table]
%                              nmfs_areas: [25×40 table]
%

% Copyright 2025 Kelly Kearney

%-------------
% Setup
%-------------

arguments
    opt.akgfmapspath {mustBeTextScalar} ="/Users/kelly.kearney/Documents/Repos/Other/akgfmaps/"
    opt.reloaddata (1,1) {mustBeNumericOrLogical} =false
    opt.dataset =[]
end

% Path to primary akgfmaps geometries

akgffol = fullfile(opt.akgfmapspath, 'inst', 'extdata');

% Path to ak_regional_mom6 supporting data (using to keep
% geopackage-derived shapefiles until I can update this function to read
% the .gpkg format natively (need to update Matlab to 2025a).

pth = fileparts(mfilename('fullpath'));
suppfol = fullfile(fileparts(pth), "supporting_data", "akgfmaps_shapefiles");

%-------------
% Read data
%-------------

akdatafile = fullfile(suppfol, 'akshapes.mat');

if opt.reloaddata

    % get_adfg_areas
    
    akmanage = readgeotable(fullfile(akgffol, 'Alaska_Marine_Management_Areas.gdb'));

    S.adfg_areas = akmanage(startsWith(akmanage.Area_Type, 'ADFG'),:);
    
    % get_base_layers
    
    S.base_layers_survey_strata = readgeotable(fullfile(suppfol, 'afsc_bts_strata_survey_strata.shp'));
    S.base_layers_survey_area   = readgeotable(fullfile(suppfol, 'afsc_bts_strata_survey_area.shp'));
    S.base_layers_inpfc_strata  = readgeotable(fullfile(suppfol, 'afsc_bts_strata_inpfc_strata.shp'));
    
    S.base_layers_longline      = readgeotable(fullfile(akgffol, 'longline_survey', 'LL_survey_Slope_and_Gullies.shp'));
    S.base_layers_ll_stations   = readgeotable(fullfile(akgffol, 'longline_survey', 'LonglineStationsActive.shp'));
    
    % get_bsierp_regions
    
    S.bsierp_regions = akmanage(startsWith(akmanage.Area_Type, 'BSIERP'),:);
    
    % get_crab_strata
    
    crab_strata = readgeotable(fullfile(suppfol,'all_crab_from_akgfmaps_grid_all_crab_from_akgfmaps_grid.shp'));
    
    S.crab_strata_RKC = crab_strata(endsWith(crab_strata.STOCK, 'RKC'),:);
    S.crab_strata_BKC = crab_strata(endsWith(crab_strata.STOCK, 'BKC'),:);
    S.crab_strata_tanner = crab_strata(startsWith(crab_strata.STOCK, 'Tanner'),:);
    S.crab_strata_snow = crab_strata(endsWith(crab_strata.STOCK, 'snow crab'),:);
    
    % ESR regions
    
    S.esr_regions = akmanage(startsWith(akmanage.Area_Type, 'Ecosystem'),:);
    
    % INPFC areas
    
    S.inpfc_goa_strata = readgeotable(fullfile(akgffol, 'goa_strata.shp'));
    S.inpfc_ai_strata = readgeotable(fullfile(akgffol, 'ai_strata.shp'));
    
    % NMFS areas
    
    S.nmfs_areas = akmanage(startsWith(akmanage.Area_Type, 'NMFS'),:);

    % Save to file

    save(akdatafile, '-struct', 'S');
end

if isempty(opt.dataset)
    S = load(akdatafile);
else
    S = load(akdatafile, opt.dataset{:});
end

end

%-------------
% Subfunctions
%-------------

function m = geo2mappolyshape(g, pcrs)
    ln = g.InternalData.VertexCoordinate2;
    lt = g.InternalData.VertexCoordinate1;
    [x,y] = projfwd(pcrs, lt, ln);
    m = mappolyshape(x,y);
    m.ProjectedCRS = pcrs;
end

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

