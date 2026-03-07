classdef cefiportalopts
%CEFIPORTALOPTS Create structure of CEFI Portal path/filename components
%
% This class holds parameters values corresponding to file/foldername
% components of the CEFI Data Portal.  
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
%   yyyymmdd:   default datestring, typically used when generating static
%               filenames
%               ["20240101"]
%
% Output variables:
%
%   C:          cefiportalopts object, structure with fields corresponding
%               short [fields] and long [fieldl] names for all path and
%               filename components, based in input options. 

    properties
        portalpath {mustBeTextScalar} ="http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/"
        region {mustBeTextScalar} ="nep"
        subdomain ="full"
        exptype {mustBeTextScalar} ="hcast"
        freq {mustBeTextScalar} ="daily"
        release {mustBeTextScalar} ="latest"
        yyyymmdd {mustBeTextScalar} ="19930101"
    end
    properties (SetAccess = private, GetAccess = public)
        regions
        regionl
        subdomains
        subdomainl
        exptypes
        exptypel
        freqs
        freql
    end

    methods
        function C = cefiportalopts(opt)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

            arguments
                opt.portalpath {mustBeTextScalar} ="http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/"
                opt.region {mustBeTextScalar} ="nep"
                opt.subdomain ="full"
                opt.exptype {mustBeTextScalar} ="hcast"
                opt.freq {mustBeTextScalar} ="daily"
                opt.release {mustBeTextScalar} ="latest"
                opt.yyyymmdd {mustBeTextScalar} ="20240101"
            end
            
            C.portalpath  = opt.portalpath;
            C.region      = opt.region;
            C.subdomain   = opt.subdomain; 
            C.exptype     = opt.exptype;
            C.freq        = opt.freq;
            C.release     = opt.release;
            C.yyyymmdd    = opt.yyyymmdd;
        end

        function C = set.region(C, val)

            tbl = [ ...
              "nwa" "northwest_atlantic"
              "nep" "northeast_pacific"
              "arc" "arctic"
              "pci" "pacific_islands"
              "glk" "great_lakes"];

            C.region = val;
            [C.regions, C.regionl] = tbllookup(tbl, val, 'region');
        end

        function C = set.subdomain(C,val)

            tbl = [...
              "full"      "full_domain"
              "aksvyreg", "ak_surveyregion_avg"];

            if ischar(val) || (isstring(val) && isscalar(val))
                if ismember(val, tbl(:))
                    C.subdomain = val;
                    [C.subdomains, C.subdomainl] = tbllookup(tbl, val, 'region');
                elseif matches(val, regexpPattern('iq\d*-\d*jq\d*-\d*'))
                    C.subdomains = val;
                    C.subdomainl = val;
                    C.subdomain  = val;
                else
                    error("subdomain not recognized");
                end
            elseif isnumeric(val) && numel(val)==4 && ...
                   allfinite(val) || all(val == floor(val), 'all')
                   
                   C.subdomains = sprintf("iq%d-%djq%d-%d", val);
                   C.subdomainl = C.subdomains;
                   C.subdomain = val;
            else
                error("subdomain not recognized");
            end

        end

        function C = set.exptype(C, val)

            tbl = [...
              "hcast"         "hindcast"
              "ss_fcast"      "seasonal_forecast"
              "ss_fcast_init" "seasonal_forecast_initialization"
              "ss_refcast"    "seasonal_reforecast"
              "dc_fcast_init" "decadal_forecast_initialization"
              "dc_fcast"      "decadal_forecast"
              "ltm_proj"      "long_term_projection"];

            C.exptype = val;
            [C.exptypes, C.exptypel] = tbllookup(tbl, val, 'exptype');

        end

        function C = set.freq(C, val)

            tbl = [...
              "day" "daily"
              "mon" "monthly"
              "ann" "yearly"
              "sta" "static"];

            C.freq = val;
            [C.freqs, C.freql] = tbllookup(tbl, val, 'freq');

        end

        function C = setopts(C, prop, val)
            arguments
                C
            end
            arguments (Repeating)
                prop {mustBeTextScalar}
                val
            end
            for ii = 1:length(prop)
                C.(prop{ii}) = val{ii};
            end
        end

        function pth = cefifolder(C, gridval)
            pth = fullfile(C.portalpath, ...
                           C.regionl, ...
                           C.subdomainl, ...
                           C.exptypel, ...
                           C.freql, ...
                           gridval, ...
                           C.release);
        end

        function fname = cefifilename(C, varname, datestr)
            fname = sprintf('%s.%s.%s.%s.%s.%s.%s.nc', ...
                        varname, ...
                        C.regions, ...
                        C.subdomains, ...
                        C.exptypes, ...
                        C.freql, ...
                        C.release, ...
                        datestr);

        end

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