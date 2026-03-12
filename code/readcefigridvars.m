function Grd = readcefigridvars(gdata, vars)
%READCEFIGRIDVARS Read grid (other) variables from one or more static files
%
% Grd = readcefigrdivars(cptopts, vars)
% Grd = readcefigrdivars(gfiles, vars)
%
% Input variables:
%
%   cpopts:     cefi portal options object corresponding to the simulation.
%               Static file names will be assumed to be the grid=raw,
%               variable=ocean_static and grid=extra, variable=ak_masks
%               files for the indicated simulation.
%
%   gfiles:     text array of explicit static file names to use
%
%   vars:       text array of variable names to read from files
%
% Output variables:
%
%   Grd:        1x1 structure array with fields corresponding to the
%               variables read from the files.

arguments
    gdata =cefiportalopts()
    vars {mustBeText} =["geolat", "geolon"]
end

% Build or verify file names

if isa(gdata, 'cefiportalopts')
    C = gdata.setopts('freq','static');
    gfiles = [...
        C.setopts('grid', 'extra').cefifilelist('ak_masks', C.yyyymmdd)
        C.setopts('grid', 'raw'  ).cefifilelist('ocean_static', C.yyyymmdd)
            ];
else
    mustBeText(gdata);
    gfiles = string(gdata);
    if ~all(cellfun(@(x) exist(x,'file'), gfiles))
        error('Specified grid files not found');
    end
end

% Read variables, iterating over variables and reading as requested.

Grd = struct;
for ii = 1:length(gfiles)
    I = ncinfo(gfiles{ii});
    vread = intersect(vars, {I.Variables.Name});
    Tmp = ncstruct(gfiles{ii}, vread{:});
    for iv = 1:length(vread)
        Grd.(vread{iv}) = Tmp.(vread{iv});
    end
    vars = setdiff(vars, vread);
end

if ~isempty(vars)
    warning('Not all requested variables found in indicated files');
end