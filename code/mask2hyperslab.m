function idx = mask2hyperslab(mask, format)
%MASK2HYPERSLAB Get hyperslab indices bounding a logical mask
%
% Input variables:
%
%   mask:   n x m logical array, size asumed to match the horizontal
%           tracer-point (i.e., h-point) model grid 
%
%   format: format to return hyperslab
%
%           "vector":   [iq1 iq2 jq1 jq2], 0-based indices of the
%                       corner-grid (q-grid) points bounding the masked
%                       region, as can be used when calling
%                       subset_from_archive.sh o extract a subregion
%
%           "ncks":     string of dimension strings that could be used to
%                       extract the relevent regions.  Note: assumes
%                       0-based indices as is the default for ncks.
%
%           "ncstruct": start-count-stride structure that can be used to
%                       extract the relevent region using the Matlab
%                       Climate Data Toolbox function ncstruct.  Note that
%                       this uses 1-based indices as is the Matlab
%                       convention.

% Copyright 2026 Kelly Kearney

arguments
    mask {mustBeNumericOrLogical, mustBeMatrix}
    format {mustBeMember(format, ["vector", "ncks", "ncstruct"])} ="vector"
end

xh = any(mask,2);
yh = any(mask,1);

% start-count-stride

switch format

    case "vector"
        idx = [find(xh,1,'first') find(xh,1,'last') ...
               find(yh,1,'first') find(yh,1,'last')] - [1 0 1 0];

    case "ncks"
        tmp = [find(xh,1,'first') find(xh,1,'last') ...
               find(yh,1,'first') find(yh,1,'last')] - 1;
        idx = sprintf("-d ih,%d,%d -d jh,%d,%d -d iq,%d,%d -d jq,%d,%d", ...
            tmp, tmp+[0 0 1 1]);

    case "ncstruct"

        idx = struct('ih', [find(xh,1) sum(xh) 1], ...
                     'jh', [find(yh,1) sum(yh) 1]);
        idx.iq = idx.ih + [0 1 0];
        idx.jq = idx.jh + [0 1 0];
end
