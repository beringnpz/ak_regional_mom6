function [hp,hpn] = plotpolyfv(T, opt)
%PLOTPOLYFV Plot face/vertex style polygons (see ak_polygons.m)
%
% Input variables:
%
%   T:  face/vertex structure (see ak_polygons Tri output)
%
% Optional variables:
%
%   val:    npoly x 1 array, where npoly is the number of polygons, i.e.,
%           max(T.fvind), in the T structure, indicating color data value
%           to be mapped to each polygon.  By default, T.color will be
%           used.
%
%   axis:   handle to axis where polygons will be plotted.  If this is a
%           map axis, the polygons will be projected appropriately.

% Copyright 2026 Kelly Kearney


    arguments
        T 
        opt.val =T.color
        opt.axis =gca
        opt.nancol (1,3) {mustBeNumeric} =nan(1,3)
    end

    if ismap(opt.axis)
        [x,y] = projfwd(getm(opt.axis), T.verts(:,2), T.verts(:,1));
        vv = [x y];
    else
        vv = T.verts;
    end
    
    hp = patch('faces', T.faces, 'vertices', vv, ...
    'facecolor', 'flat', 'facealpha', 1.0, 'edgecolor', 'none', ...
    'cdata', opt.val(T.fvind));

    if ~all(isnan(opt.nancol))
        mask = isnan(opt.val(T.fvind));
        hpn = patch('faces', T.faces(mask,:), 'vertices', vv, ...
            'facecolor', opt.nancol, 'facealpha', 1.0, 'edgecolor', 'none');

    end
        


end