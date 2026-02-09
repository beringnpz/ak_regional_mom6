function [Idx, regnum] = readindexdata(simname, opt)
%READINDEXDATA Reads data from survey-region-average files
%
% [Idx, regnum] = readindexdata(simname, opt)
%
% This function reads in data from the survey-region-average files, holding
% spatially-averaged timeseries from the hindcast, plus persistence
% forecast and climatology for that hindcast.
%
% Input variables:
%
%   simname:    name of simulation, used to locate output data files.  The
%               index data will be read from files under
%               <datafol>/<simname>/Level3/surveyregionavg/
%
% Optional input variables (passed as parameter/value pairs, default in []):
%
%   datafol:    CEFI data folder path.  Default is the path returned by the
%               cefidatafol.m function
%
%   vars:       variable names for which index data will be returned.  Can
%               be any of the variables in the file, or one of the
%               following additional ones that will be calculated based on
%               those: pH, omega.  If empty, all variables found in file
%               will be returned. [''] 
%
%   regnum:     region numbers for which index data will be returned.  If
%               NaN, all regions found in file will be used [NaN] 
%
% Output variables:
%
%   Tbl:        1 x 1 structure with the following fields:
%
%               Hc:     timetable of hindcast index data, with variables
%                       (columns) corresponding to the requested variables
%                       and each holding a nt x nreg array of values. 
%
%               Fc:     timetable of persistence forecast data, with
%                       similar structure as Hc
%
%               Clim:   timetable of climatological hindcast index data,
%                       with similar structure as Hc
%
%   regnum:     nreg x 1 array, region number values corresponding to data
%               in Tbl (if input regnum was NaN or not provided, this will
%               correspond to all regions found in the survey region files;
%               otherwise, it will match the input value).

% Copyright 2025 Kelly Kearney

arguments
    simname {mustBeTextScalar}
    opt.datafol {mustBeTextScalar} =cefidatafolpath
    opt.vars {mustBeText} =''
    opt.regnum {mustBeNumeric} =NaN
end

% Survey-region-based index data files

fnamehc = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*selected_daily*svyreg.nc'));
fnamefc = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*forecast*svyreg.nc'));
fnamecm = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*daily_clim*svyreg.nc'));

% Check variables in file vs requested ones

ftmp = fullfile(fnamehc(1).folder, fnamehc(1).name);
I = ncinfo(ftmp);
isv = arrayfun(@(X) isequal({X.Dimensions.Name}, {'time', 'surveyregion'}), I.Variables);

varsinfile = {I.Variables(isv).Name};

if isempty(opt.vars)
    opt.vars = varsinfile;
end
opt.vars = string(opt.vars);

% Adjust reading for potential extra variables (pH, omega)

varsread = setdiff(opt.vars, ["pH", "omega"]);
if ismember("pH", opt.vars)
    varsread = union(varsread, "btm_htotal");
end
if ismember("omega", opt.vars)
    varsread = union(varsread, ["btm_co3_ion","btm_co3_sol_arag"]);
end

% vars = ["tob", "tos", "btm_htotal", "btm_o2", "cpool0p0", "cpool2p0", "btm_co3_ion", "btm_co3_sol_arag"];
% vlong = ["Bottom temp. (\circC)", "SST (\circC)", "Bottom pH", "Bottom O_2 (mmol/kg)"];

% Check regions

regnum = ncread(ftmp, 'mask_survey_area');
if isscalar(opt.regnum) && isnan(opt.regnum)
    opt.regnum = regnum;
end

[tf,loc] = ismember(opt.regnum, regnum);
if ~all(tf)
    error('Unrecognized region number');
end

% Read data

Idx.Hc   = readandconcat(fnamehc, varsread);
Idx.Fc   = readandconcat(fnamefc, varsread);
Idx.Clim = readandconcat(fnamecm, varsread);

% A few adjustments and additions

if ismember("pH", opt.vars)
    Idx.Hc.pH = -log10(Idx.Hc.btm_htotal./1.25);
    Idx.Fc.pH = -log10(Idx.Fc.btm_htotal./1.25);
    Idx.Clim.pH = -log10(Idx.Clim.btm_htotal./1.25);
end

if ismember("omega", opt.vars)
    Idx.Hc.omega = Idx.Hc.btm_co3_ion./Idx.Hc.btm_co3_sol_arag;
    Idx.Fc.omega = Idx.Fc.btm_co3_ion./Idx.Fc.btm_co3_sol_arag;
    Idx.Clim.omega = Idx.Clim.btm_co3_ion./Idx.Clim.btm_co3_sol_arag;
end

% Convert to timetable and extract only specified regions

Idx = structfun(@(X) idxstruct2table(X, loc, opt.vars), Idx, 'uni', 0);

end

% Subfunction: read selected data into a structure

function A = readandconcat(fname, vars)

    fname = fullfile({fname.folder}, {fname.name});

    Tmp = arrayfun(@(x) ncstruct(x, vars{:}, 'time'), fname);
    for iv = 1:length(vars)
        A.(vars{iv}) = cat(1, Tmp.(vars{iv}));
    end
    A.time = cat(1, Tmp.time);
    A.t = ncdateread(fname{1}, 'time', A.time);

end

% Subfunction: convert index structure to timetable

function B = idxstruct2table(A, loc, vplt)
    B = table2timetable(struct2table(A));
    B = removevars(B, 'time');
    B = B(:,vplt);
    B = varfun(@(x) x(:,loc), B);
    B.Properties.VariableNames = strrep(B.Properties.VariableNames, 'Fun_', '');
end


