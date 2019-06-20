function output = importtxt(path_to_file, delimiter, startRow)
%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NCwb19\Overview.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2018/09/12 20:06:39

%% Initialize variables.
%path_to_file = 'C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NCwb19\Overview.txt';
%delimiter = ',';
%startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(path_to_file,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.



%% Create output variable
output = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

