function a = geotable2struct(b, opt)
%GEOTABLE2STRUCT Convert geospatial table to a geographic data structure
%
% a = geotable2struct(b)
%
% Input variables:
%
%   b:  table containing, at minimum, a Shape column with either
%       mappolyshape or geopolyshape objects
%
% Optional input variables:
%
%   proj:   if true, reverse-project mappolyshape object coordinates to
%           lat/lon based on each shape's CRS [false]
%
%   simp:   if non-NaN, apply line simplification to shapes using the
%           specified tolerance (in units of degrees for geopolyshapes or
%           the CRS unit for mappolyshapes
%
%   wrap:   180: wrap any longitude coordinates to [-180 180]
%           360: wrap any longitude coordinates to [0 360]
%           0:   leave as is
%           [0]
%
%   poly:   if true, create a polyshape object [false]
%
% Output variables:
%
%   a:  n x 1 geographic data structure, where n is the height of table b.  
%       If b held geopolyshapes, a will contain Lat/Lon fields with
%       coordinate data; if b help mappolyshapes, it will contain x/y data.

arguments
    b 
    opt.proj (1,1) {mustBeNumericOrLogical} =false
    opt.simp (1,1) =NaN
    opt.wrap (1,1) {mustBeNumeric, mustBeMember(opt.wrap, [0 180 360])} =0
    opt.poly (1,1) {mustBeNumericOrLogical} =false
end

if ~isgeotable(b)
   error('First input must be a geotable')
end

a = table2struct(b);

for ii = 1:length(a)
    x = a(ii).Shape.InternalData.VertexCoordinate1;
    y = a(ii).Shape.InternalData.VertexCoordinate2;
    idx = a(ii).Shape.InternalData.IndexOfLastVertex;

    n = idx - [1 idx(1:end-1)+1]+1;

    x = mat2cell(x, 1, n);
    y = mat2cell(y, 1, n);
    x = cellfun(@(z) [z NaN], x, 'uni', 0);
    y = cellfun(@(z) [z NaN], y, 'uni', 0);
    x = cat(2, x{:});
    y = cat(2, y{:});

    if isa(a(ii).Shape, 'geopolyshape')
        a(ii).Geometry = 'Polygon';
        a(ii).Lat = x;
        a(ii).Lon = y;
    elseif isa(a(ii).Shape, 'mappolyshape')
        a(ii).Geometry = 'Polygon';
        
        if ~isnan(opt.simp)
            xysimp = dpsimplify([x; y]', opt.simp);
            a(ii).Y = [xysimp(:,2); NaN]';
            a(ii).X = [xysimp(:,1); NaN]';
        else
            a(ii).Y = y;
            a(ii).X = x;
        end

        if opt.proj
            [a(ii).Lat, a(ii).Lon] = projinv(a(ii).Shape.ProjectedCRS, x, y);
        end
    end
    if opt.wrap == 360
        a(ii).Lon = wrapTo360(a(ii).Lon);
    elseif opt.wrap == 180
        a(ii).Lon = wrapTo180(a(ii).Lon);
    end
    if opt.poly
        a(ii).poly = polyshape(a(ii).Lon, a(ii).Lat);
    end
end

