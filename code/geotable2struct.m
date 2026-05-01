function a = geotable2struct(b)
%VERTEXDATA Convert geospatial table to a geographic data structure

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
            a(ii).Lat = x;
            a(ii).Lon = wrapTo360(y);
        elseif isa(a(ii).Shape, 'mappolyshape')
            a(ii).Y = y;
            a(ii).X = x;
        end
    end

end