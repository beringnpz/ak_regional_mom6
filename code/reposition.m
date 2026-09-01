function reposition(hobj, scale, shift, anchor)
%REPOSITION Scale and/or shift a graphics object's position
%
% Input variables:
%
%   hobj:   handle to graphics object
%   scale:  1x2 array, scale factor to apply to horizontal and vertical
%           positions, respectively
%   shift:  1x2 array: shift to apply to anchor point position, in units
%           used by hobj Position attribute
%   anchor: anchor point indicating corner or side to anchor in place
%           during rescaling  

    arguments
        hobj (1,1)
        scale (1,2) {mustBeNumeric}
        shift (1,2) {mustBeNumeric}
        anchor {mustBeTextScalar, mustBeMember(anchor, {'nw','n','ne','e','se','s','sw','w'})}
    end

    if ~isgraphics(hobj)
        error('First input (hobj) must be valid graphics object');
    end

    pos = get(hobj, 'Position');
    
    w = pos(3)*scale(1);
    h = pos(4)*scale(2);
    
    A = table(...
        [...
            pos(1)          pos(2)+pos(4)
            pos(1)+pos(3)/2 pos(2)+pos(4)
            pos(1)+pos(3)   pos(2)+pos(4)
            pos(1)+pos(3)   pos(2)+pos(4)/2
            pos(1)+pos(3)   pos(2)
            pos(1)+pos(3)/2 pos(2)
            pos(1)          pos(2)
            pos(1)          pos(2)+pos(4)/2], ...
        [...
            0       -h
            -w/2    -h
            -w      -h
            -w      -h/2
            -w      0
            -w/2    0
            0       0
            0       -h/2], ...
        'rownames', {'nw','n','ne','e','se','s','sw','w'}, ...
        'variablenames', {'refxy', 'shift'});
    
    corner = A{anchor,'refxy'} + shift + A{anchor,'shift'};
    newpos = [corner w h];
    
    set(hobj, 'position', newpos);

end