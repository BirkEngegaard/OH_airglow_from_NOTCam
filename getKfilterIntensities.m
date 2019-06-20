function [q97counts, relerr97, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = getKfilterIntensities(file, masterDark, masterFlat,showFig)
q97counts = -1;
relerr97 = -1;
specificDark = -1;
exposureTime = -1;
obsTime = -1;
airmass = -1;
SAT = -1;
rampMode = -1;

%% Read in data
wd = 'D:\NOTArchive\'; % working directory
day = file(1:6);
path = [wd, file(1:6), '\']; % Where the raw data is

%% Make spec
data = fitsread(strcat(path,file,'.fits'),'Image');
S=fitsinfo(strcat(path,file,'.fits'));
k=S.PrimaryData.Keywords;
[objindex,~] = find(strcmp(k,'OBJECT'));
[r, ~] = size(objindex); % rows, columns
if r == 2
    objindex = objindex(1);
end
if contains(k{objindex,2}, 'Ha_flat')
    return
end

[obsindex,~] = find(strcmp(k,'OBS_MODE'));
if ~strcmp(k{obsindex,2}, 'Spectroscopy')
    return
end
[aindex1, ~] = find(strcmp(k,'APERTUR'));
[aindex2, ~] = find(strcmp(k,'NCAPRNM'));
if ~isempty(aindex1)
    if ~strcmp(k{aindex1,2}, '128mu WF slit')
        return;
    end
elseif ~isempty(aindex2)
    if ~strcmp(k{aindex2,2}, '128mu WF slit')
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
    warning('Couldn''t determine filter, return -1')
    return;
end

%% Get exposure time and choose the observation data point to be in the middle of observation
expMode = S.PrimaryData.Keywords{find(strcmp(k,'EXPMODE')), 2}; % Exposure mode e.g. 'frames 3.6 3'
dexpMode = ['d', expMode];
expMode = split(expMode); % split string into array for easier handling

try
    exposureTime = str2double(expMode{2})*str2double(expMode{3}); % Total time = time of ramp * num of ramps
catch
    warning('Couldn''t determine exposure time. Default value of 10.8 used')
    exposureTime = 10.8;
end

%figpath = sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s-%s_figs\\',filter, file); % To save figures for debug/study
figpath = sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\data\\%s-%s\\',filter, file); % To save figures for debug/study
if not ( exist(figpath))
    mkdir(figpath)
end


[q97counts, relerr97, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = prepare_K(path,file,filter,figpath, S, data, masterDark, masterFlat, showFig);

end