function datafol = cefidatafolpath
%CEFIDATAFOLPATH Returns the path of the CEFI data folder
%
% datafol = cefidatafolpath
%
% The ak_regional_mom6 code is designed such that Alaska-region datasets
% are stored separately from the code repository.  This function can be
% used to set and retrieve the path name.
%
% If the ak_regional_mom6/simulation_data/data_folder.txt file does not
% exist, calling this function will present an interactive folder dialog
% box where the user can choose the location of the data folder; the
% data_folder.txt file will be updated with the chosen path.  If the file
% does exist, the function simply return the path saved within that file.  

% Copyright 2025 Kelly Kearney

datafolfile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'simulation_data', 'data_folder.txt');
if ~exist(datafolfile, 'file')
    datafol = uigetdir(fileparts(datafolfile), 'Select ak_regional_mom6-related CEFI Data folder');
    fid = fopen(datafolfile, 'wt');
    fprintf(fid, '%s', datafol);
    fclose(fid);
else
    datafol = fileread(datafolfile);
end