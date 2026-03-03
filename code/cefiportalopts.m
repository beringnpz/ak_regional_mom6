function C = cefiportalopts(opt)
%CEFIPORTALOPTS Create structure of CEFI Portal path/filename components
%
% Input variables (passed as parameter/value pairs, default in []):
%
%   portalpath: base path for CEFI Portal (or copy) location.  By default,
%               this will be the THREDDS Data Server for the web-based CEFI
%               Data Portal.  
%               ["http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/"]
%
%   region:     short or long region name:
%                 "nwa" "northwest"
%                 "nep" "northeast_pacific"
%                 "arc" "arctic"
%                 "pci" "pacific_islands"
%                 "glk" "great_lakes"
%               ["nep"]
%
%   subdomain:  short or long subdomain name: 
%                 "full" "full_domain" "iq#-#jq#-#"
%               or [iqmin iqmax jqmin jqmax] vector of corner-point indices
%               describing a subdomain
%               ["full"]
%
%   exptype:    short or long experiment type name
%                 "hcast"         "hindcast"
%                 "ss_fcast"      "seasonal_forecast"
%                 "ss_fcast_init" "seasonal_forecast_initialization"
%                 "ss_refcast"    "seasonal_reforecast"
%                 "dc_fcast_init" "decadal_forecast_initialization"
%                 "dc_fcast"      "decadal_forecast"
%                 "ltm_proj"      "long_term_projection"
%               ["hcast"]
%
%   freq:       short or long frequency name:
%                 "day" "daily"
%                 "mon" "monthly"
%                 "ann" "yearly"
%                 "sta" "static"
%               ["daily"]
%
%   release:    release name
%
% Output variables:
%
%   C:          1x1 structure with fields corresponding short [fields] and
%               long [fieldl] names for all path and filename components,
%               based in input options. 

% Copyright 2026 Kelly Kearney

arguments
    opt.portalpath {mustBeTextScalar} ="http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/"
    opt.region {mustBeTextScalar} ="nep"
    opt.subdomain ="full"
    opt.exptype {mustBeTextScalar} ="hcast"
    opt.freq {mustBeTextScalar} ="daily"
    opt.release {mustBeTextScalar} ="latest"
end

T.region =[ ...
  "nwa" "northwest"
  "nep" "northeast_pacific"
  "arc" "arctic"
  "pci" "pacific_islands"
  "glk" "great_lakes"];

T.subdomain = [...
  "full" "full_domain"];

T.exptype = [...
  "hcast"         "hindcast"
  "ss_fcast"      "seasonal_forecast"
  "ss_fcast_init" "seasonal_forecast_initialization"
  "ss_refcast"    "seasonal_reforecast"
  "dc_fcast_init" "decadal_forecast_initialization"
  "dc_fcast"      "decadal_forecast"
  "ltm_proj"      "long_term_projection"
    ];

T.freq = [...
  "day" "daily"
  "mon" "monthly"
  "ann" "yearly"
  "sta" "static"];

% Simple matches

C.portalpath = opt.portalpath;
C.release = opt.release;

% Subdomain

try mustBeTextScalar(opt.subdomain)
    if ismember(opt.subdomain, T.subdomain(:))
        % will be handled by long/short lookup
    elseif matches(opt.subdomain, regexpPattern('iq\d*-\d*jq\d*-\d*'))
        C.subdomains = opt.subdomain;
        C.subdomainl = opt.subdomain;
        T = rmfield(T, "subdomain");
    else
        error("subdomain not recognized");
    end
catch
    if isnumeric(opt.subdomain) && numel(opt.subdomain)==4 && ...
       allfinite(opt.subdomain) || all(opt.subdomain == floor(opt.subdomain), 'all')
       C.subdomains = sprintf("iq%d-%djq%d-%d", opt.subdomain);
       C.subdomainl = C.subdomains;
       T = rmfield(T, "subdomain");
    else
        error("subdomain not recognized");
    end
end

% Long/short lookup

for vv = string(fieldnames(T))'
    [C.(vv+"s"), C.(vv+"l")] =  tbllookup(T.(vv), opt.(vv), vv);
end

end

function [sname, lname] = tbllookup(tbl, x, v)
%TBLLOOKUP 2-column string table lookup

    [tf,loc] = ismember(x, tbl);
    if tf
        [idx,~] = ind2sub(size(tbl),loc);
        sname = tbl(idx,1);
        lname = tbl(idx,2);
    else
        error("Option (%s) not recognized; %s must match one of: %s", x, v, join(tbl(:), ", "));
    end

end