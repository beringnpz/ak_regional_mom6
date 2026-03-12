function [Idx, regnum] = readindexdata(cpopts, opt)
%READINDEXDATA Reads data from survey-region-average files
%
% [Idx, regnum] = readindexdata(opt)
%
% This function reads in data from the survey-region-average files, holding
% spatially-averaged timeseries from the hindcast, plus persistence
% forecast and climatology for that hindcast.
%
% Optional input variables (passed as parameter/value pairs, default in []):
%
%   cpopts:     cefi portal options object corresponding to the simulation
%
%   vars:       variable names for which index data will be returned.  Can
%               be any of the variables in the files, or one of the
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
    cpopts (1,1) {mustBeA(cpopts, "cefiportalopts")} =cefiportalopts()
    opt.vars {mustBeText} =''
    opt.regnum {mustBeNumeric} =NaN
end

% Check variables in survey-region-based index data files

favail = cpopts.setopts('grid', 'extra').cefifilelist('*','*');
[~, fshort] = fileparts(favail);

vmain = unique(regexp(fshort, '^[^\.]*', 'match', 'once'));

varsinfile = cell(size(vmain));
for iv = 1:length(vmain)
    I = ncinfo(favail{iv});
    isv = arrayfun(@(X) ~isempty(X.Dimensions) && isequal({X.Dimensions.Name}, {'time', 'surveyregion'}), I.Variables);
    varsinfile{iv} = {I.Variables(isv).Name};
end
varsinfile = cat(2, varsinfile{:});

if isempty(opt.vars)
    opt.vars = varsinfile;
end
opt.vars = string(opt.vars);

% Adjust reading for potential extra variables (pH, omega) or tagalong
% variables (cpoolX.X)

varsread = setdiff(opt.vars, ["pH", "omega"]);
if ismember("pH", opt.vars)
    varsread = union(varsread, "btm_htotal");
end
if ismember("omega", opt.vars)
    varsread = union(varsread, ["btm_co3_ion","btm_co3_sol_arag"]);
end

% Check regions

regnum = ncread(favail{1}, 'mask_survey_area');
if isscalar(opt.regnum) && isnan(opt.regnum)
    opt.regnum = regnum;
end

[tf,loc] = ismember(opt.regnum, regnum);
if ~all(tf)
    error('Unrecognized region number');
end

% Read data

Idx.Hc   = struct;
Idx.Fc   = struct;
Idx.Clim = struct;

for iv = 1:length(varsread)

    if startsWith(varsread{iv}, 'cpool')
        vfile = "tob";
    else
        vfile = varsread{iv};
    end
   
    flisthc = cpopts.setopts('grid','extra').cefifilelist(vfile, '*');
    flistfc = cpopts.setopts('grid','extra').cefifilelist("fcpersist_"+vfile, '*');
    flistcm = cpopts.setopts('grid','extra').cefifilelist("clim_"+vfile, '*');

    Idx.Hc   = readandconcat(Idx.Hc,   flisthc, varsread{iv}, iv==1);
    Idx.Fc   = readandconcat(Idx.Fc,   flistfc, varsread{iv}, iv==1);
    Idx.Clim = readandconcat(Idx.Clim, flistcm, varsread{iv}, iv==1);

end

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

function A = readandconcat(A, fname, vv, tflag)

    if tflag
        Tmp = arrayfun(@(x) ncstruct(x, vv, 'time'), fname);
    else
        Tmp = arrayfun(@(x) ncstruct(x, vv), fname);
    end

    if tflag
        A.time = cat(1, Tmp.time);
        A.t = ncdateread(fname{1}, 'time', A.time);
    end

    A.(vv) = cat(1, Tmp.(vv));

end

% Subfunction: convert index structure to timetable

function B = idxstruct2table(A, loc, vplt)
    B = table2timetable(struct2table(A));
    B = removevars(B, 'time');
    B = B(:,vplt);
    B = varfun(@(x) x(:,loc), B);
    B.Properties.VariableNames = strrep(B.Properties.VariableNames, 'Fun_', '');
end


