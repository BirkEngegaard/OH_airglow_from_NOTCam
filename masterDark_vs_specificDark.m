%% Analysis of intensity variation using master dark and specific dark
% Runtime = 25 minutes for K band
%% load variables
OHintensityData = load('OHintensity_data_HK_220319.mat');

H_specDarks = OHintensityData.specificDarks;
K_specDarks = OHintensityData.K_specificDarks;

H_used_files = OHintensityData.used_files;
K_used_files = OHintensityData.K_used_files;

intensity31 = OHintensityData.intensity31;
intensity42 = OHintensityData.intensity42;
intensity53 = OHintensityData.intensity53;
intensity64 = OHintensityData.intensity64;
intensity97 = OHintensityData.intensity97;

%% Master darks and flats
dpath = 'C:\Users\birke\Documents\master\scripting\masterDarks\';
dark = dlmread(strcat(dpath,'masterDark_raw_median.txt')); % dark in ADU/sec

%% H filter
hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));

%% K filter
K_hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
K_flat = dlmread(strcat(hpath, 'K_masterFlat_raw_median.txt'));

%% Get names of files with specific darks
H_names = H_used_files(H_specDarks == 1);
K_names = K_used_files(K_specDarks == 1);

%% Get intensities using specific darks
intensity31_specDark = intensity31(H_specDarks == 1);
intensity42_specDark = intensity42(H_specDarks == 1);
intensity53_specDark = intensity53(H_specDarks == 1);
intensity64_specDark = intensity64(H_specDarks == 1);
intensity97_specDark = intensity97(K_specDarks == 1);

%% Initialize variables for intensities using master dark
intensity31_masterDark = zeros(1,length(H_names));
intensity42_masterDark = zeros(1,length(H_names));
intensity53_masterDark = zeros(1,length(H_names));
intensity64_masterDark = zeros(1,length(H_names));
intensity97_masterDark = zeros(1,length(K_names));

intensity31_masterDark_relerr = zeros(1,length(H_names));
intensity42_masterDark_relerr = zeros(1,length(H_names));
intensity53_masterDark_relerr = zeros(1,length(H_names));
intensity64_masterDark_relerr = zeros(1,length(H_names));
intensity97_masterDark_relerr = zeros(1,length(K_names));

%% Scan through H band with master dark
progress = 0;
total_progress = length(H_names) + length(K_names);
w = waitbar(progress/total_progress, 'Scanning H band');

if true
for i = 1:length(H_names)
    [Int31, relerr31,Int42,relerr42,Int53,relerr53,Int64,relerr64, SAT,...
        specificDark, expTime, obsTime, airmass]...
        = SinglePrepare_dev(H_names{i}, dark, flat, false);
    
    intensity31_masterDark(i) = Int31;
    intensity31_masterDark_relerr(i) = relerr31;
    intensity42_masterDark(i) = Int42;
    intensity42_masterDark_relerr(i) = relerr42;
    intensity53_masterDark(i) = Int53;
    intensity53_masterDark_relerr(i) = relerr53;
    intensity64_masterDark(i) = Int64;
    intensity64_masterDark_relerr(i) = relerr64;
    
    progress = progress + 1;
    waitbar(progress/total_progress,w,'Scanning H band')
end
end
%% Scan through K band with master dark

for i = 1:length(K_names)
[int97, relerr97, specificDark, expTime, obsTime, airmass]...
    = getKfilterIntensities(K_names{i}, dark, K_flat,false);
    
    intensity97_masterDark(i) = int97;
    intensity97_masterDark_relerr(i) = relerr97;
    
    progress = progress + 1;
    waitbar(progress/total_progress,w,'Scanning K band')
end

close(w)

difference31 = intensity31_masterDark - intensity31_specDark;
rel_difference31 = difference31./intensity31_specDark;
difference42 = intensity42_masterDark - intensity42_specDark;
rel_difference42 = difference42./intensity42_specDark;
difference53 = intensity53_masterDark - intensity53_specDark;
rel_difference53 = difference53./intensity53_specDark;
difference64 = intensity64_masterDark - intensity64_specDark;
rel_difference64 = difference64./intensity64_specDark;

difference97 = intensity97_masterDark - intensity97_specDark;
rel_difference97 = difference97./intensity97_specDark;

H_rel_diff = [rel_difference31 rel_difference42 rel_difference53 rel_difference64];
tot_rel_diff = [rel_difference97 H_rel_diff];

fitdist(transpose(rel_difference97), 'normal')
figure;
histogram(rel_difference97)
title('relative error 97')

figure;
plot(rel_difference97)
