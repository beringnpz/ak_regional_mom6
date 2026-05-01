function h = plotmom6akmap(A, X, opt)
%PLOTMOM6AKMAP Plots a map of a horizontal slice of MOM6 output
%
% h = plotmom6akmap(A, X)
% h = plotmom6akmap(A, X, param, value, ...)
%
% This function provides a quick plot of MOM6 spatial data, tailored to
% Alaska regional uses.
%
% Input variables:
%
%   A:  structure of plotting details returned from akmapprep function
%
%   X:  structure of spatial slice data to be used for plotted maps, as
%       returned by readmom6mapslice function
%
% Optional input variables (passed as parameter/value pairs):
%
%   var:        variable name, as found in the X data structure (and used
%               to look up visualization data via the ak_variable_info
%               table
%
%   vartype:    'anomaly' or 'value', used to set colormaps appropriately
%
%   mapprops:   cell array of parameter/value pairs to pass to boxworldmap
%
%   cbloc:      colorbar position.  If empty, no colorbar will be added.
%
%   addsyrveypoly:  logical scalar, true to plot outline of the Eastern
%               Bering Sea groundfish survey

% Copyright 2026 Kelly Kearney

arguments
    A (1,1) % akmapprep structure
    X (1,1) % readmom6mapslice structure
    opt.var {mustBeTextScalar} ='tob'
    opt.vartype {mustBeTextScalar} ='value'
    opt.mapprops ={} % to pass to boxworldmap 
    opt.cbloc ='west'
    opt.addsurveypoly =false
end

h.bm = A.addboxmap(opt.mapprops);
h.pc = pcolorpad(A.xc, A.yc, padend(X.(opt.var)));
shading flat;
uistack(h.pc, 'bottom');

% Add land borders 

h.bor = plot(A.bx, A.by, 'color', rgb('gray'));

if opt.addsurveypoly
    h.svy = plot(A.sx, A.sy, 'k');
end

% Set colormap and add colorbar

V = ak_variable_info;

switch opt.vartype
    case 'anomaly'
        set(gca, 'clim', V{opt.var, 'limanom'}, 'colormap',cmocean('balance'), 'layer', 'top');
    case 'value'
        set(gca, 'clim', V{opt.var, 'limmap'}, 'colormap', V{opt.var, 'cmap'}{1}, 'layer', 'top');
end
if ~isempty(opt.cbloc)
    h.cb = colorbar(opt.cbloc);
end

end

function x = padend(x)
%PADEND Add trailing row and column of NaNs to a 2D array

    x = [x nan(size(x,1),1); nan(1,size(x,2)+1)];

end