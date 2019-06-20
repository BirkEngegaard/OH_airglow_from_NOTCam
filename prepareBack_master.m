function [dark,hot,flat,bad, useMasterDark, nightMaster] = prepareBack_master(dpath, file,filter, masterDark, masterFlat)

[useMasterDark, dlist, dexposureTime, nightMaster] = findDarks(file);
%% Remove this after masterDark_vs_specificDark
% useMasterDark = 1;
% nightmaster = 0;
%%
if ~nightMaster % Only do this specific darks
    [N, ~] = size(dlist);
    if N > 1
        dlist(1,:) = [];
    elseif N > 2
        dlist(1,:) = [];
        dlist(2,:) = [];
    elseif N > 3
        dlist(1,:) = [];
        dlist(2,:) = [];
        dlist(3,:) = [];
    end
end
%% darks
if useMasterDark
    if strcmp(masterFlat, 'none')
        dpath = 'C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\masterDarks\';
        dark = dlmread(strcat(dpath,'masterDark_raw_median.txt')); % dark in ADU/sec
    else
        dark = masterDark;
    end
elseif ~nightMaster  
try
    dark = fitsread(strcat(dpath,dlist(1,:),'.fits'),'Image');
catch
    errstr = sprintf('Darkfields error for file: %s',file);
    error(errstr)
end
% figure;
% contourf(dark,'LineColor','none');
% pcolor(dark); shading interp; 
% colormap('hot(600)')
% colorbar;
% title("Dark field 1, unprocessed")

for i=2:size(dlist)
     dark(:,:,i) = fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image');
end

% if strcmp(file, 'NCyi280117') || strcmp(file, 'NCyi280118')...
%         || strcmp(file, 'NCyi280119') || strcmp(file, 'NCyi280120')...
%         || strcmp(file, 'NCyi280121')
%     for i=2:size(dlist)
%         dark(:,:,i) = 0.04*fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image'); % scale down darkfield
%     end
% elseif strcmp(file, 'NCyi280131') || strcmp(file, 'NCyi280132')...
%         || strcmp(file, 'NCyi280133') || strcmp(file, 'NCyi280134')...
%         || strcmp(file, 'NCyi280135')
%     for i=2:size(dlist)
%         dark(:,:,i) = 0.1*fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image'); % scale down darkfield
%     end
% else
%     for i=2:size(dlist)
%         dark(:,:,i) = fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image');
%     end
% end

%%
% dS = fitsinfo(strcat(dpath,dlist(1,:),'.fits'));
% dk = dS.PrimaryData.Keywords;
% expMode = dS.PrimaryData.Keywords{find(strcmp(dk,'EXPMODE')), 2}; % Exposure mode e.g. 'frames 3.6 3'
% dexpMode = ['d', expMode];
% expMode = split(expMode, ' '); % split string into array for easier handling
% try
%     dexposureTime = str2double(expMode{2})*str2double(expMode{3}); % Total time = time of ramp * num of ramps
% catch
%     warning('Couldn''t determine exposure time. Default value of 10.8 used')
%     dexposureTime = 10.8;
% end
%%
   
%dS = fitsinfo(strcat(dpath,dlist(2,:),'.fits'));
%dk = dS.PrimaryData.Keywords;
%2 combine with median for each pixel
dark = median(dark,3);
%dark = mean(dark, 3);
dark = dark/dexposureTime;
elseif nightMaster
    for i=1:size(dlist)
        dark(:,:,i) = fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image')/dexposureTime(i);
    end
    dark = median(dark,3);
end
%%

% if strcmp(file, 'NCsa100199') || strcmp(file, 'NCsa100188') || strcmp(file, 'NCsa100189') || strcmp(file, 'NCsa100190')...
%         || strcmp(file, 'NCsa100191') || strcmp(file, 'NCsa100196') || strcmp(file, 'NCsa100197') || strcmp(file, 'NCsa100198')
%     dpath = 'D:\NOTArchive\NCsa10\';
%     dexposureTime = 16;
%     dlist = ['NCsa100435';];
%     for i=1:size(dlist)
%         dark(:,:,i) = fitsread(strcat(dpath,dlist(i,:),'.fits'),'Image');
%     end
%     dark = median(dark,3);
%     dark = dark/dexposureTime;
% end

bad = fitsread('bad_zero_sci.fits');
bad(bad > 0)=1;
dark(bad == 1) = 10001;
%3 smooth over 10x10 pixels (median) -> dark
hot = [handleBackground(dark(1:512,1:512)),handleBackground(dark(1:512,513:1024));handleBackground(dark(513:1024,1:512)),handleBackground(dark(513:1024,513:1024))];
hot = [handleBackground2(hot(1:512,1:512)),handleBackground2(hot(1:512,513:1024));handleBackground2(hot(513:1024,1:512)),handleBackground2(hot(513:1024,513:1024))];

%dark = dark - hot;
dark(:,512)=0;
dark(bad == 1) = 0;

%% Hot 0 if over 6 sigma, else 1 for each pixel -> hot
hot(hot>550/10.8)=0;
%dark(dark<-10000)=0;
hotaverage = mean(mean(hot));
hotsigma = std2(hot);
hot(abs(hot)>abs(hotaverage)+6*hotsigma)=-100001;%0.06
hot(hot>=-100000)=1;
hot(hot<-100000)=0;
%disp(sum(sum(abs(hot),1))/1024/1024);

%% flatfield:
if strcmp(filter, 'H')
    if strcmp(masterFlat, 'none')
        hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
        %flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));
        flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));
    else
        flat = masterFlat;
    end
elseif strcmp(filter, 'K')
    if strcmp(masterFlat, 'none')
        hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
        %flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));
        flat = dlmread(strcat(hpath, 'K_masterFlat_raw_median.txt'));
    else
        flat = masterFlat;
    end
elseif strcmp(filter, 'J')
    if strcmp(masterFlat, 'none')
        hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
        %flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));
        flat = dlmread(strcat(hpath, 'J_masterFlat_raw_median.txt'));
    else
        flat = masterFlat;
    end
end
%hpath = 'C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\masterFlats\';
%combine heliums with median
%flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2009_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2010_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2011_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2012_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2013_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2014_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2015_raw.txt'));
%flat = dlmread(strcat(hpath, 'meanFlat2016_raw.txt'));

%substract zero -> put to high number
flat = flat - dark;
flat(bad == 1) = 0;%1000001;
%substract hot -> put to high number
flat(hot == 0) = 0;%1000001;
%normalize that the average pixel is 1 in the middle (corners are gonna explode)
%flat = flat/mean(mean(flat)); % Lets try not doing this
flat(flat==0)=1;

end