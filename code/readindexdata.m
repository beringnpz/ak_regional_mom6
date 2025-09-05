function [Idx, regnum] = readindexdata(simname, opt)

arguments
    simname {mustBeTextScalar}
    opt.staticname {mustBeTextScalar} =simname+"_ocean_static_ak.nc"
    opt.datafol {mustBeTextScalar} =cefidatafolpath
end

% Target variables

vars = ["tob", "tos", "btm_htotal", "btm_o2", "cpool0p0", "cpool2p0", "btm_co3_ion", "btm_co3_sol_arag"];
% vlong = ["Bottom temp. (\circC)", "SST (\circC)", "Bottom pH", "Bottom O_2 (mmol/kg)"];

% Survey-based index data

fnamehc = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*selected_daily*svyreg.nc'));
fnamefc = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*forecast*svyreg.nc'));
fnamecm = dir(fullfile(opt.datafol, simname, 'Level3', 'surveyregionavg', '*daily_clim*svyreg.nc'));

[Idx.Hc, regnum] = readandconcat(fnamehc, vars);
Idx.Fc = readandconcat(fnamefc, vars);
Idx.Clim = readandconcat(fnamecm, vars);

% A few adjustments and additions

% Hc.pH = -log10(Hc.btm_htotal./1.25);
% Fc.pH = -log10(Fc.btm_htotal./1.25);
% 
% Hc.btm_o2 = Hc.btm_o2*1e6;
% Fc.btm_o2 = Fc.btm_o2*1e6;
% 
% Hc.omega = Hc.btm_co3_ion./Hc.btm_co3_sol_arag;
% Fc.omega = Fc.btm_co3_ion./Fc.btm_co3_sol_arag;



end

function [A, regnum] = readandconcat(fname, vars)

    fname = fullfile({fname.folder}, {fname.name});

    Tmp = arrayfun(@(x) ncstruct(x, vars{:}, 'time'), fname);
    for iv = 1:length(vars)
        A.(vars{iv}) = cat(1, Tmp.(vars{iv}));
    end
    A.time = cat(1, Tmp.time);
    A.t = ncdateread(fname{1}, 'time', A.time);

    if nargin > 1
        regnum = ncread(fname{1}, 'mask_survey_area');
    end

end

