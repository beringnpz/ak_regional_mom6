function archivepng(h, folder, filename, varargin)
%ARCHIVEPNG Save copy of figure to a png file
%
% This quick wrapper function calls export_fig to save a figure to a png
% image file.  It only creates an image if the file does not already exist.
% This is primarily intended so I can quickly save higher-resolution
% versions of figures to an external archive folder, apart from the
% snapshots that are added to the quarto book automatically.
%
% archivepng(h, folder, filename)
% archivepng(h, folder, filename, '-p1', '-p2', ...)
%
% Input variables:
%
%   h:          handle to figure
%
%   folder:     path to folder where image will be saved
%
%   filename:   filename, without extension
%
%   -p1, ...:   optional inputs to pass to export_fig   

% Copyright 2026 Kelly Kearney

fname = fullfile(folder, filename);

if ~exist([fname '.png'], 'file')
    export_fig(fname, h, '-png', varargin{:});
end