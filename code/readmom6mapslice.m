function X = readmom6mapslice(C, ttarget, opt)
%READMOM6MAPSLICE Reads single time step of 2D MOM6 variable
%
% Input variables:
%
%   C:          cefi portal options object corresponding to the simulation
%
%   ttarget:    scalar datetime, date to extract from output
%
%   vars:       string array of variable names to extract from the dataset.
%               Currently, these must be 2D variables.
%
%   vartype:    type of data to read:
%               'value':    hindcast value
%               'anomaly':  anomaly from climatology
%               'fcpersis': persistence forecast value

% Copyright 2026 Kelly Kearney

arguments
    C (1,1) {mustBeA(C, "cefiportalopts")}
    ttarget (1,1) {mustBeA(ttarget, "datetime")}
    opt.vars {mustBeText} ="tob"
    opt.vartype {mustBeTextScalar} ='value'
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

for iv = 1:length(varsread)

    switch opt.vartype
        case 'anomaly'
            fname = C.setopts('freq','daily','grid','extra').cefifilelist("anom_"+varsread{iv}, year(ttarget)+"*");
        case 'value'
            fname = C.setopts('freq','daily','grid','raw').cefifilelist(varsread{iv}, year(ttarget)+"*");
        case 'fcpersist'
            fname = C.setopts('freq','daily','grid','extra').cefifilelist("fcpersist_"+varsread{iv}, year(ttarget)+"*");
        otherwise
            error('Unrecognized variable type (%s); should be value or anomaly', opt.vartype)
    end

    t = ncdateread(fname, 'time');
    [dt,imin] = min(abs(ttarget - t));
    if dt > days(5)
        warning('Possible gap: %s is closest time found', t(imin));
    end
    Tmp = ncstruct(fname, varsread{iv}, struct('time', [imin 1 1]));
    X.(varsread{iv}) = Tmp.(varsread{iv});
end

if ismember("pH", opt.vars) % TODO: doesn't work for anomaly
    X.pH = -log10(X.btm_htotal./1.25);
end

if ismember("omega", opt.vars) % TODO: doesn't work for anomaly
    X.omega = X.btm_co3_ion./X.btm_co3_sol_arag;
end

if ismember("btm_o2", opt.vars)
    X.btm_o2 = X.btm_o2*1e6;
end

