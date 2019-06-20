function [q31counts, relerr31, q42counts, relerr42,q53counts,relerr53,q64counts,relerr64, SAT,...
    specificDark, exposureTime, obsTime, airmass, rampMode]...
    = SinglePrepare_dev(file, masterDark, masterFlat,showFig)
specificDark = -1;
exposureTime = -1;
obsTime = -1;
airmass = -1;
SAT = -1;
rampMode = -1;
%close all;
%% Read in data
wd = 'D:\NOTArchive\'; % working directory
%wd = 'C:\Users\birke\Documents\master\scripting\';


%file = 'NCwb191099'; % image
%file = 'NCwb191101'; % K
%file = 'NCwb191102';
%file = 'NCwb191103';
%file = 'NCwb191104';
%file = 'NCwb191105'; % H
%file = 'NCwb191106';
%file = 'NCwb191107'; % Interesting dump
%file = 'NCwb191108';
%file = 'NCwb191109'; % J filter

%file = 'NCza230395';

day = file(1:6);
path = [wd, file(1:6), '\\']; % Where the raw data is
dpath = path; %Where is the dark image (if not in the same folder)
hpath = path; %Where is the flatfield image (if not in the same folder)
trans = 42; %which transition

% read overview file

%overview = importtxt([wd, file(1:6),'\Overview.txt'], ',', 2); % Imports Overview.txt
 
%filter1 = 'K'; %which filter, for debug
%filter = findFilter(overview, file); % Scans Overview.txt and finds filter used for file
%if (filter1 ~= filter)
%    error('Filter name error!')
%end


%% Make spec
data = fitsread(strcat(path,file,'.fits'),'Image');
S=fitsinfo(strcat(path,file,'.fits'));
k=S.PrimaryData.Keywords;
[obsindex,~] = find(strcmp(k,'OBS_MODE'));
if ~strcmp(k{obsindex,2}, 'Spectroscopy')
    q31counts = -1;
    q42counts = -1;
    q53counts = -1;
    q64counts = -1;
    SAT = -1;
    relerr31 = -1;
    relerr42 = -1;
    relerr53 = -1;
    relerr64 = -1;
    return
end

[aindex1, ~] = find(strcmp(k,'APERTUR'));
[aindex2, ~] = find(strcmp(k,'NCAPRNM'));
if ~isempty(aindex1)
    if ~strcmp(k{aindex1,2}, '128mu WF slit')
        q31counts = -2;
        q42counts = -2;
        q53counts = -2;
        q64counts = -2;
        SAT = -2;
        relerr31 = -2;
        relerr42 = -2;
        relerr53 = -2;
        relerr64 = -2;
        return;
    end
elseif ~isempty(aindex2)
    if ~strcmp(k{aindex2,2}, '128mu WF slit')
        q31counts = -2;
        q42counts = -2;
        q53counts = -2;
        q64counts = -2;
        SAT = -2;
        relerr31 = -2;
        relerr42 = -2;
        relerr53 = -2;
        relerr64 = -2;
        return;
    end
end

[findex,~] = find(strcmp(k,'NCFLTNM2'));
[findex2,~] = find(strcmp(k,'FILT2'));
if ~isempty(findex)
    filter = k{findex,2};
elseif ~isempty(findex2)
    filter = k{findex2,2};
else
    q31counts = -1;
    q42counts = -1;
    q53counts = -1;
    q64counts = -1;
    SAT = -1;
    relerr31 = -1;
    relerr42 = -1;
    relerr53 = -1;
    relerr64 = -1;
    warning('Couldn''t determine filter, return -1')
    return;
end



%% Get exposure time and choose the observation data point to be in the middle of observation
expMode = S.PrimaryData.Keywords{find(strcmp(k,'EXPMODE')), 2}; % Exposure mode e.g. 'frames 3.6 3'
dexpMode = ['d', expMode];
expMode = split(expMode, ' '); % split string into array for easier handling

try
    exposureTime = str2double(expMode{2})*str2double(expMode{3}); % Total time = time of ramp * num of ramps
catch
    warning('Couldn''t determine exposure time. Default value of 10.8 used')
    exposureTime = 10.8;
end

%[flatfields, darkframes, Hfilter, master] = masterfileanalysis(filter, dexpMode);

hlist = [];
dlist = [];
% for i = 1:size(flatfields)
%     if strcmp(flatfields{i}(1:6), day)
%         hlist = [hlist; flatfields{i}];
%     end
% end
% 
% for i = 1:size(darkframes)
%     if strcmp(darkframes{i}(1:6), day)
%         dlist = [dlist; darkframes{i}];
%     end
% end
% 
% [m, n] = size(dlist);
% if m > 1
% dlist(1,:) = [];
% end
% 
% if strcmp(file, 'NCyi280117') || strcmp(file, 'NCyi280118')...
%         || strcmp(file, 'NCyi280119') || strcmp(file, 'NCyi280120')...
%         || strcmp(file, 'NCyi280121') || strcmp(file, 'NCyi280119')...
%         || strcmp(file, 'NCyi280131') || strcmp(file, 'NCyi280132')...
%         || strcmp(file, 'NCyi280133') || strcmp(file, 'NCyi280134')...
%         || strcmp(file, 'NCyi280135')
%     dlist = ['NCyi280002'; 'NCyi280003'; 'NCyi280004'; 'NCyi280005'; 'NCyi280006'];
% end
    


%hlist = findFlatfields(path, filter);
%dlist = findDarkfields(path);

%figpath = sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s-%s_figs\\',filter, file); % To save figures for debug/study
figpath = sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\data\\%s-%s\\',filter, file); % To save figures for debug/study
if not ( exist(figpath))
    mkdir(figpath)
end

%hlist = ['NCwb191089';'NCwb191090']; %J
hlist = ['NCwb191093';'NCwb191094']; %H
%hlist = ['NCwb191097';'NCwb191098']; %K
dlist = ['NCwb191174';'NCwb191175';'NCwb191176';'NCwb191177'];

[q31counts,relerr31,q42counts,relerr42,q53counts,relerr53,q64counts,relerr64,SAT,specificDark, exposureTime, obsTime, airmass, rampMode]...
    = prepare(path,file,dpath,dlist,hpath,hlist,trans,filter,figpath, S, data, masterDark, masterFlat, showFig);

plotting(path, figpath, file, filter, dlist, hlist, trans)

end