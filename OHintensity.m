% Runtime = 2 hours

%% Which "modules" to run
runHfilter = false;
runKfilter = false;
runJfilter = true;

%% Load file names
H_files = load('C:\Users\birke\Documents\master\scripting\Hfilterdata_correct.mat');
H_files = H_files.Hfilterdata;

K_files = load('C:\Users\birke\Documents\master\scripting\Kfilterdata.mat');
K_files = K_files.Kfilterdata;

J_files = load('C:\Users\birke\Documents\master\scripting\Jfilterdata.mat');
J_files = J_files.Jfilterdata;

%% H and K data
%OHintensityData = load('OHintensity_data_HK_070319.mat');

%% Latex as default text interpreter for figures
set(0,'DefaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%% Master darks and flats
dpath = 'C:\Users\birke\Documents\master\scripting\masterDarks\';
dark = dlmread(strcat(dpath,'masterDark_raw_median.txt')); % dark in ADU/sec

%% H filter
hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
flat = dlmread(strcat(hpath, 'masterFlat_raw_median.txt'));

%% K filter
K_hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
K_flat = dlmread(strcat(hpath, 'K_masterFlat_raw_median.txt'));

%% J filter
J_hpath = 'C:\Users\birke\Documents\master\scripting\masterFlats\';
J_flat = dlmread(strcat(hpath, 'J_masterFlat_raw_median.txt'));

%%
%Hdays = load('C:\Users\birke\Documents\master\scripting\Hfilter_days.mat');
%Hdays = Hdays.Hfilter_days; % Column cell array
%day = 'NCwb19';
%day = 'NCyi28';
%day = 'NCre23';
%day = 'NCsk28';
%day = 'NCyi27';
%day = 'NCvf12';
%day = 'NCwf22';
%day = 'NCyg03';
%day = 'NCyi04';
%day = 'NCyi28';
%day = 'NCre22';
%day = 'NCwf22';
%day = 'NCyg03';
%day = 'NCti30';
%day = 'NCub12';
% day = 'NCyg03';
%H_files = findHfiles(day);
%clear day;

%% H filter init
[r,c] = size(H_files);
N = r;

obstime = cell(1,N);
intensity31 = zeros(1,N);
relative_error31 = zeros(1,N);
intensity42 = zeros(1,N);
relative_error42 = zeros(1,N);
intensity53 = zeros(1,N);
relative_error53 = zeros(1,N);
intensity64 = zeros(1,N);
relative_error64 = zeros(1,N);
specificDarks = zeros(1,N);
exposureTime = zeros(1,N);
used_files = cell(1,N);
unused_files = cell(1,N);
time = cell(1,N);
time(1,:) = {0};
used_files(1,:) = {0};
unused_files(1,:) = {0};
obstime(1,:) = {0};
nightOfYear = zeros(1,N);
hPostSunSet = zeros(1,N);
H_sdep = zeros(1,N);
H_airmass = zeros(1,N);
H_LST = cell(1,N);
H_LST(1,:) = {0};
H_rampExp = zeros(1,N); % ramp exposure time
H_rampNum = zeros(1,N); % Number of ramps
k = 0;
j = 0;

%% K filter init
[K_r,K_c] = size(K_files);
K_N = K_r;

K_obstime = cell(1,K_N);
intensity97 = zeros(1,K_N);
relative_error97 = zeros(1,K_N);
K_specificDarks = zeros(1,K_N);
K_exposureTime = zeros(1,K_N);
K_used_files = cell(1,K_N);
K_unused_files = cell(1,K_N);
K_time = cell(1,K_N);
K_time(1,:) = {0};
K_used_files(1,:) = {0};
K_unused_files(1,:) = {0};
K_obstime(1,:) = {0};
K_nightOfYear = zeros(1,K_N); 
K_hPostSunSet = zeros(1,K_N); % Hours post sunset at NOT
K_sdep = zeros(1,K_N); % Solar depression angle
K_airmass = zeros(1,K_N);
K_LST = cell(1,K_N);
K_LST(1,:) = {0};
K_rampExp = zeros(1,K_N); % ramp exposure time
K_rampNum = zeros(1,K_N); % Number of ramps
K_k = 0;
K_j = 0;

%% J filter init
[J_r,J_c] = size(J_files);
J_N = J_r;

J_obstime = cell(1,J_N);
intensity74 = zeros(1,J_N);
relative_error74 = zeros(1,J_N);
intensity85 = zeros(1,J_N);
relative_error85 = zeros(1,J_N);
J_specificDarks = zeros(1,J_N);
J_exposureTime = zeros(1,J_N);
J_used_files = cell(1,J_N);
J_unused_files = cell(1,J_N);
J_time = cell(1,J_N);
J_time(1,:) = {0};
J_used_files(1,:) = {0};
J_unused_files(1,:) = {0};
J_obstime(1,:) = {0};
J_nightOfYear = zeros(1,J_N); 
J_hPostSunSet = zeros(1,J_N); % Hours post sunset at NOT
J_sdep = zeros(1,J_N); % Solar depression angle
J_airmass = zeros(1,J_N);
J_LST = cell(1,J_N);
J_LST(1,:) = {0};
J_rampExp = zeros(1,J_N); % ramp exposure time
J_rampNum = zeros(1,J_N); % Number of ramps
J_k = 0;
J_j = 0;

%% For data spread evaluation
evaluateHdates = false;
if evaluateHdates
    for i = 1:N
        [data, S] = getFitsData(H_files{i});
        [oindex, ~] = find(strcmp(S.PrimaryData.Keywords,'DATE-OBS'));
        obsTime = datetime(replace(S.PrimaryData.Keywords{oindex, 2}, 'T', ' '), 'Format', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
        obstime{i} = obsTime;
    end
    histogram2([day([obstime{:}], 'dayofyear')], [year([obstime{:}])],[365,8], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year')
    ylabel('Year')
end

%% Spectral analysis
%% H Filter

if runHfilter
wb = waitbar(0,'(1/2) Scanning H filter files');
for i = 1:N %(length(H_files)-N):length(H_files)
    [Int31, relerr31,Int42,relerr42,Int53,relerr53,Int64,relerr64,...
        SAT, specificDark, expTime, obsTime, airmass, rampMode]...
        = SinglePrepare_dev(H_files{i}, dark, flat, false);
    %[Int, relerr, SAT] = SinglePrepare_dev(H_files(i,:), dark, flat);
    if (Int31 ~= -1) && (Int31 ~= -2)...
            && ~strcmp(H_files{i}, 'NCyi040817')... % This one has ridiculous errorbars
            && ~strcmp(H_files{i},'NCwf220793')... % This one is negative, probably bad
            && ~strcmp(H_files{i}, 'NCyi040786')... % Also negative, probably bad. Note NCyi04 being thrown twice
            && ~strcmp(H_files{i},  'NCyi040819')... % Same as above
            && ~strcmp(H_files{i}, 'NCsf090430') % Very noisy, throw it away

        k = k + 1;
        intensity31(k) = Int31;
        intensity42(k) = Int42;
        intensity53(k) = Int53;
        intensity64(k) = Int64;
        
        if specificDark ~= 1
            darkError = 0.1; % 
        else
            darkError = 0;
        end

        relerr31 = abs(relerr31) + darkError; % +10% Because of darkfields;
        relerr42 = abs(relerr42) + darkError;
        relerr53 = abs(relerr53) + darkError;
        relerr64 = abs(relerr64) + darkError;
        relative_error31(k) = abs(relerr31);
        relative_error42(k) = abs(relerr42);
        relative_error53(k) = abs(relerr53);
        relative_error64(k) = abs(relerr64);
        
        %time{k} = datetime(SAT,'Format','yyyy/MM/dd HH:mm:ss');
        time{k} = obsTime; % UTC
        H_LST{k} = SAT; % Local solar time = Solar apparent time
        used_files{k} = H_files{i};
        specificDarks(k) = specificDark; % 1 = specific dark, 0 = master dark, 2 = master dark from the same night
        exposureTime(k) = expTime;
        nightOfYear(k) = day(datetime(sprintf('%d/%d/%s',H_files{i}(3) + 1894, H_files{i}(4) - 96, H_files{i}(5:6)), 'InputFormat', 'yyyy/MM/dd'), 'dayofyear');
        obsYear = H_files{i}(3) + 1894;
        obsMonth = H_files{i}(4) - 96;
        obsDay = str2num(H_files{i}(5:6));
        nightDate = datetime(obsYear, obsMonth, obsDay);
        [~,SSET,~] = sunrise(28.75,-17.883,0,0, datenum(nightDate));
        SSET = datetime(SSET, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
        hPostSunSet_k = hours(datetime(obsTime, 'TimeZone', 'UTC') - SSET);
        hPostSunSet(k) = hPostSunSet_k;
        [~, El] = SolarAzEl(obsTime,28.75,-17.883,3);
        H_sdep(k) = -El; % degrees
        H_airmass(k) = airmass;
        H_rampExp(k) = rampMode(1);
        H_rampNum(k) = rampMode(2);
        
    else
        j = j + 1;
        unused_files{j} = H_files{i};
    end
    if i/N < 0.5
        waitbar(i/N,wb,'(1/2) Scanning H filter files');
    elseif i/N > 0.5 && i/N < 0.8
        waitbar(i/N,wb,'(1/2) Scanning H filter files, the worst is over!');
    elseif i/N > 0.8 && i/N <0.9
        waitbar(i/N,wb,'(1/2) Scanning H filter files, just a few more now...')
    elseif i/N > 0.9
        waitbar(i/N,wb, '(1/2) Scanning H filter files, don''t quit on us now!')
    end
end
close(wb)
end

%% K filter spectral data gathering

if runKfilter
    K_wb = waitbar(0, '(2/2) Scanning K filter files');
    for i = 1:K_N
        [int97, relerr97, specificDark, expTime, obsTime, airmass, SAT, rampMode] = getKfilterIntensities(K_files{i}, dark, K_flat,false);
        if (int97 ~= -1) && (int97 ~= -2)...
                && ~strcmp(K_files{i}, 'NCqk160348')...% This one is bad 
                && ~strcmp(K_files{i}, 'NCqk160349')...% This one is very noisy
                && ~strcmp(K_files{i}, 'NCqk160350')...% Same as above
                && ~strcmp(K_files{i}, 'NCqk160351')...% Same as above
                && ~strcmp(K_files{i}, 'NCvi060814')... %Can't see any lines in raw data despite 60s exptime
                && ~strcmp(K_files{i}, 'NCvj210295')... % Too noisy
                && ~strcmp(K_files{i}, 'NCvj210296')...
                && ~strcmp(K_files{i}, 'NCvj210297')...
                && ~strcmp(K_files{i}, 'NCvj210298')...
                && ~strcmp(K_files{i}, 'NCwf220694')...
                && ~strcmp(K_files{i}, 'NCwf220697')...
                && ~strcmp(K_files{i}, 'NCwf220698')...
                && ~strcmp(K_files{i}, 'NCwf220701')...
                && ~strcmp(K_files{i}, 'NCwf220702')...
                && ~strcmp(K_files{i}, 'NCxh090177')...
                && ~strcmp(K_files{i}, 'NCxh090178')...
                && ~strcmp(K_files{i}, 'NCxh090179')...
                && ~strcmp(K_files{i}, 'NCxh090180')...
                && ~strcmp(K_files{i}, 'NCyg030314')...
                && ~strcmp(K_files{i}, 'NCyg030315')...
                && ~strcmp(K_files{i}, 'NCyg030316')...
                && ~strcmp(K_files{i}, 'NCyg030317')...
                && ~strcmp(K_files{i}, 'NCyg030364')...
                && ~strcmp(K_files{i}, 'NCyg030365')...
                && ~strcmp(K_files{i}, 'NCyg030366')...
                && ~strcmp(K_files{i}, 'NCyg030239')...
                && ~strcmp(K_files{i}, 'NCyg030240')...
                && ~strcmp(K_files{i}, 'NCyg030241')...
                && ~strcmp(K_files{i}, 'NCyg030242')...
                && ~strcmp(K_files{i}, 'NCyg040744')...
                && ~strcmp(K_files{i}, 'NCyg040745')...
                && ~strcmp(K_files{i}, 'NCyg040746')...
                && ~strcmp(K_files{i}, 'NCyg040747')...
                && ~strcmp(K_files{i}, 'NCyi020239')...
                && ~strcmp(K_files{i}, 'NCyi020240')...
                && ~strcmp(K_files{i}, 'NCyi020241')...
                && ~strcmp(K_files{i}, 'NCyi020242')...
                && ~strcmp(K_files{i}, 'NCyi040744')...
                && ~strcmp(K_files{i}, 'NCyi040745')...
                && ~strcmp(K_files{i}, 'NCyi040746')...
                && ~strcmp(K_files{i}, 'NCyi040747')...
                && ~strcmp(K_files{i}, 'NCyi040813')...
                && ~strcmp(K_files{i}, 'NCyi040814')...
                && ~strcmp(K_files{i}, 'NCyi040815')...
                && ~strcmp(K_files{i}, 'NCyi040816')...
                && ~strcmp(K_files{i}, 'NCyi270370')...
                && ~strcmp(K_files{i}, 'NCyi280190')...
                && ~strcmp(K_files{i}, 'NCyi280191')...
                && ~strcmp(K_files{i}, 'NCyi280192')...
                && ~strcmp(K_files{i}, 'NCyi280193')...
                && ~strcmp(K_files{i}, 'NCyi280194')
            K_k = K_k + 1;
            intensity97(K_k) = int97;
            
            if specificDark ~= 1
                darkError = 0.1; % 
            else
                darkError = 0;
            end
            relerr97 = abs(relerr97) + darkError; % +10% Because of darkfields;
            relative_error97(K_k) = abs(relerr97);
            K_time{K_k} = obsTime; % UTC
            K_LST{K_k} = SAT; % Local solar time, solar apparent time
            K_used_files{K_k} = K_files{i};
            K_specificDarks(K_k) = specificDark;
            K_exposureTime(K_k) = expTime;
            K_nightOfYear(K_k) = day(datetime(sprintf('%d/%d/%s',K_files{i}(3) + 1894, K_files{i}(4) - 96, K_files{i}(5:6)), 'InputFormat', 'yyyy/MM/dd'), 'dayofyear');
            obsYear = K_files{i}(3) + 1894;
            obsMonth = K_files{i}(4) - 96;
            obsDay = str2num(K_files{i}(5:6));
            nightDate = datetime(obsYear, obsMonth, obsDay);
            [~,SSET,~] = sunrise(28.75,-17.883,0,0, datenum(nightDate));
            SSET = datetime(SSET, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
            hPostSunSet_k = hours(datetime(obsTime, 'TimeZone', 'UTC') - SSET);
            K_hPostSunSet(K_k) = hPostSunSet_k;
            [~, El] = SolarAzEl(obsTime,28.75,-17.883,3);
            K_sdep(K_k) = -El; % degrees
            K_airmass(K_k) = airmass;
            K_rampExp(K_k) = rampMode(1);
            K_rampNum(K_k) = rampMode(2);

        else
            K_j = K_j + 1;
            K_unused_files{K_j} = K_files{i};
        end
        
        if i/K_N < 0.5
            waitbar(i/K_N,K_wb,'(2/2) Scanning K filter files');
        elseif i/K_N > 0.5 && i/N < 0.8
            waitbar(i/K_N,K_wb,'(2/2) Scanning K filter files, the worst is over!');
        elseif i/K_N > 0.8 && i/N <0.9
            waitbar(i/K_N,K_wb,'(2/2) Scanning K filter files, just a few more now...')
        elseif i/K_N > 0.9
            waitbar(i/K_N,K_wb, '(2/2) Scanning K filter files, don''t quit on us now!')
        end
    end
    close(K_wb)
end

%% J filter spectral data gathering
J_wb = waitbar(0, '(3/3) Scanning J filter files');
if runJfilter
    for i = 1:J_N
        [int74, relerr74, int85, relerr85, specificDark, expTime, obsTime, airmass, SAT, rampMode] = getJfilterIntensities(J_files{i}, dark, J_flat,false);
        if (int74 ~= -1) && (int74 ~= -2)
            J_k = J_k + 1;
            intensity74(J_k) = int74;
            intensity85(J_k) = int85;
            
            if specificDark ~= 1
                darkError = 0.1; % 
            else
                darkError = 0;
            end
            relerr74 = abs(relerr74) + darkError; % +10% Because of darkfields;
            relative_error74(J_k) = abs(relerr74);
            relerr85 = abs(relerr85) + darkError; % +10% Because of darkfields;
            relative_error85(J_k) = abs(relerr85);
            J_time{J_k} = obsTime; % UTC
            J_LST{J_k} = SAT; % Local solar time, solar apparent time
            J_used_files{J_k} = J_files{i};
            J_specificDarks(J_k) = specificDark;
            J_exposureTime(J_k) = expTime;
            J_nightOfYear(J_k) = day(datetime(sprintf('%d/%d/%s',J_files{i}(3) + 1894, J_files{i}(4) - 96, J_files{i}(5:6)), 'InputFormat', 'yyyy/MM/dd'), 'dayofyear');
            obsYear = J_files{i}(3) + 1894;
            obsMonth = J_files{i}(4) - 96;
            obsDay = str2num(J_files{i}(5:6));
            nightDate = datetime(obsYear, obsMonth, obsDay);
            [~,SSET,~] = sunrise(28.75,-17.883,0,0, datenum(nightDate));
            SSET = datetime(SSET, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
            hPostSunSet_k = hours(datetime(obsTime, 'TimeZone', 'UTC') - SSET);
            J_hPostSunSet(J_k) = hPostSunSet_k;
            [~, El] = SolarAzEl(obsTime,28.75,-17.883,3);
            J_sdep(J_k) = -El; % degrees
            J_airmass(J_k) = airmass;
            J_rampExp(J_k) = rampMode(1);
            J_rampNum(J_k) = rampMode(2);

        else
            J_j = J_j + 1;
            J_unused_files{J_j} = J_files{i};
        end
        
        if i/J_N < 0.5
            waitbar(i/J_N,J_wb,'(3/3) Scanning J filter files');
        elseif i/J_N > 0.5 && i/N < 0.8
            waitbar(i/J_N,J_wb,'(3/3) Scanning J filter files, the worst is over!');
        elseif i/J_N > 0.8 && i/N <0.9
            waitbar(i/J_N,J_wb,'(3/3) Scanning J filter files, just a few more now...')
        elseif i/J_N > 0.9
            waitbar(i/J_N,J_wb, '(3/3) Scanning J filter files, don''t quit on us now!')
        end
    end
    close(J_wb)
end

%% H filter | cutaway zero-tail
intensity31 = intensity31(1:k);
relative_error31 = relative_error31(1:k);
intensity42 = intensity42(1:k);
relative_error42 = relative_error42(1:k);
intensity53 = intensity53(1:k);
relative_error53 = relative_error53(1:k);
intensity64 = intensity64(1:k);
relative_error64 = relative_error64(1:k);
used_files = used_files(1:k);
unused_files = unused_files(1:k);
time = time(1:k);
H_time = time;
H_LST = H_LST(1:k);
specificDarks = transpose(specificDarks(1:k));
exposureTime = exposureTime(1:k);
H_exposureTime = exposureTime(1:k);
hPostSunSet = hPostSunSet(1:k);
H_hpss = hPostSunSet;
nightOfYear = nightOfYear(1:k);
H_nightOfYear = nightOfYear;
H_noy = H_nightOfYear;
H_sdep = H_sdep(1:k);
H_airmass = H_airmass(1:k);
H_rampExp = H_rampExp(1:k);
H_rampNum = H_rampNum(1:k);

%% K filter | cutaway zero-tail
intensity97 = intensity97(1:K_k);
relative_error97 = relative_error97(1:K_k);
K_used_files = K_used_files(1:K_k);
K_unused_files = K_unused_files(1:K_j);
K_time = K_time(1:K_k);
K_LST = K_LST(1:K_k);
K_specificDarks = transpose(K_specificDarks(1:K_k));
K_exposureTime = K_exposureTime(1:K_k);
K_hPostSunSet = K_hPostSunSet(1:K_k);
K_hpss = K_hPostSunSet;
K_nightOfYear = K_nightOfYear(1:K_k);
K_noy = K_nightOfYear;
K_sdep = K_sdep(1:K_k);
K_airmass = K_airmass(1:K_k);
K_rampExp = K_rampExp(1:K_k);
K_rampNum = K_rampNum(1:K_k);

%% J filter | cutaway zero-tail
intensity74 = intensity74(1:J_k);
relative_error74 = relative_error74(1:J_k);
intensity85 = intensity85(1:J_k);
relative_error85 = relative_error85(1:J_k);
J_used_files = J_used_files(1:J_k);
J_unused_files = J_unused_files(1:J_j);
J_time = J_time(1:J_k);
J_LST = J_LST(1:J_k);
J_specificDarks = transpose(J_specificDarks(1:J_k));
J_exposureTime = J_exposureTime(1:J_k);
J_hPostSunSet = J_hPostSunSet(1:J_k);
J_hpss = J_hPostSunSet;
J_nightOfYear = J_nightOfYear(1:J_k);
J_noy = J_nightOfYear;
J_sdep = J_sdep(1:J_k);
J_airmass = J_airmass(1:J_k);
J_rampExp = J_rampExp(1:J_k);
J_rampNum = J_rampNum(1:J_k);

disp('Stopping script in line 470)')
return

%% H filter | Absolute error calcs
abserr31 = abs(relative_error31.*intensity31);
abserr42 = abs(relative_error42.*intensity42);
abserr53 = abs(relative_error53.*intensity53);
abserr64 = abs(relative_error64.*intensity64);
time_vec = transpose([time{:}]);

%% K filter | Absolute error calc
abserr97 = abs(relative_error97.*intensity97);
K_time_vec = transpose([K_time{:}]);



%% H filter | Relative intensities, with errors, and weights
i42 = intensity42(:)./intensity31(:);
relerr_i42 = sqrt(relative_error31(:).^2 + relative_error42(:).^2);
w42 = 1./(relerr_i42(:).^2);
i53 = intensity53(:)./intensity31(:);
relerr_i53 = sqrt(relative_error31(:).^2 + relative_error53(:).^2);
w53 = 1./(relerr_i53(:).^2);
i64 = intensity64(:)./intensity31(:);
relerr_i64 = sqrt(relative_error31(:).^2 + relative_error64(:).^2);
w64 = 1./(relerr_i64(:).^2);

%% H filter | Intensity difference relative to 31, with errors, and weights
i42_diff = intensity42(:) - intensity31(:);
relerr_i42_diff = (1./i42_diff).*sqrt(abserr42(:).^2 + abserr31.^2);
w42_diff = 1./(relerr_i42_diff(:).^2);
i53_diff = intensity53(:) - intensity31(:);
relerr_i53_diff = (1./i53_diff).*sqrt(abserr53(:).^2 + abserr31.^2);
w53_diff = 1./(relerr_i53_diff(:).^2);
i64_diff = intensity64(:) - intensity31(:);
relerr_i64_diff = (1./i64_diff).*sqrt(abserr64(:).^2 + abserr31.^2);
w64_diff = 1./(relerr_i64_diff(:).^2);

%% Yearly variation
%% H filter | annual variation
intensity31_year = zeros(1,12);
err31_year = zeros(1,12);
intensity31_year_relerr30 = zeros(1,12);
err31_year_relerr30 = zeros(1,12);
intensity31_year_weighted = zeros(1,12);
err31_year_weighted = zeros(1,12);
intensity42_year = zeros(1,12);
intensity42_year_relerr30 = zeros(1,12);
intensity42_year_weighted = zeros(1,12);
err42_year = zeros(1,12);
err42_year_relerr30 = zeros(1,12);
err42_year_weighted = zeros(1,12);
intensity53_year = zeros(1,12);
intensity53_year_relerr30 = zeros(1,12);
intensity53_year_weighted = zeros(1,12);
err53_year = zeros(1,12);
err53_year_relerr30 = zeros(1,12);
err53_year_weighted = zeros(1,12);
intensity64_year = zeros(1,12);
intensity64_year_relerr30 = zeros(1,12);
intensity64_year_weighted = zeros(1,12);
err64_year = zeros(1,12);
err64_year_relerr30 = zeros(1,12);
err64_year_weighted = zeros(1,12);

if runHfilter
for i = 1:12
    %% Global mean
    intensity31_year(i) = mean(intensity31(month([time{:}]) == i));
    err31_year(i) = 1/length(abserr31(month([time{:}]) == i)).*sqrt((sum([abserr31(month([time{:}]) == i).^2])));
    intensity42_year(i) = mean(intensity42(month([time{:}]) == i));
    err42_year(i) = 1/length(abserr42(month([time{:}]) == i)).*sqrt((sum([abserr42(month([time{:}]) == i).^2])));
    intensity53_year(i) = mean(intensity53(month([time{:}]) == i));
    err53_year(i) = 1/length(abserr53(month([time{:}]) == i)).*sqrt((sum([abserr53(month([time{:}]) == i).^2])));
    intensity64_year(i) = mean(intensity64(month([time{:}]) == i));
    err64_year(i) = 1/length(abserr64(month([time{:}]) == i)).*sqrt((sum([abserr64(month([time{:}]) == i).^2])));

    %% weighted mean
    intensity31_year_weighted(i) = sum([intensity31(month([time{:}]) == i).*(1./(relative_error31(month([time{:}]) == i).^2))])...
        /sum(1./relative_error31(month([time{:}]) == i).^2);
    
    err31_year_weighted(i) = 1/sum([1./relative_error31(month([time{:}]) == i).^2])*...
        sqrt((sum([((1./relative_error31(month([time{:}]) == i).^2)).*abserr31(month([time{:}]) == i)])));
    
    intensity42_year_weighted(i) = sum([intensity42(month([time{:}]) == i).*(1./(relative_error42(month([time{:}]) == i).^2))])...
        /sum(1./relative_error42(month([time{:}]) == i).^2);
    err42_year_weighted(i) = 1/sum([1./relative_error42(month([time{:}]) == i).^2])*...
        sqrt((sum([(1./relative_error42(month([time{:}]) == i).^2).*abserr42(month([time{:}]) == i)])));
    intensity53_year_weighted(i) = sum([intensity53(month([time{:}]) == i).*(1./(relative_error53(month([time{:}]) == i).^2))])...
        /sum(1./relative_error53(month([time{:}]) == i).^2);
    err53_year_weighted(i) = 1/sum([1./relative_error53(month([time{:}]) == i).^2])*...
        sqrt((sum([(1./relative_error53(month([time{:}]) == i).^2).*abserr53(month([time{:}]) == i)])));
    intensity64_year_weighted(i) = sum([intensity64(month([time{:}]) == i).*(1./(relative_error64(month([time{:}]) == i).^2))])...
        /sum(1./relative_error64(month([time{:}]) == i).^2);
    err64_year_weighted(i) = 1/sum([1./relative_error64(month([time{:}]) == i).^2])*...
        sqrt((sum([(1./relative_error64(month([time{:}]) == i).^2).*abserr64(month([time{:}]) == i)])));    
    
end
end

%% K filter 
%% K filter | Annual variation
intensity97_year = zeros(1,12);
err97_year = zeros(1,12);
intensity97_year_weighted = zeros(1,12);
err97_year_weighted = zeros(1,12);

if runKfilter
    for i = 1:12
        %% Global mean
        intensity97_year(i) = mean(intensity97(month([K_time{:}]) == i));
        err97_year(i) = 1/length(abserr97(month([K_time{:}]) == i)).*sqrt((sum([abserr97(month([K_time{:}]) == i).^2])));
        
        %% Weighted mean
        intensity97_year_weighted(i) = sum([intensity97(month([K_time{:}]) == i).*(1./(relative_error97(month([K_time{:}]) == i).^2))])...
        /sum(1./relative_error97(month([K_time{:}]) == i).^2);
    
        err97_year_weighted(i) = 1/sum([1./relative_error97(month([K_time{:}]) == i).^2])*...
        sqrt((sum([((1./relative_error97(month([K_time{:}]) == i).^2)).*abserr97(month([K_time{:}]) == i)])));
    
        err97_year_weighted2(i) = 1/sum([1./relative_error97(month([K_time{:}]) == i).^2])*...
            sqrt(sum((intensity97(month([K_time{:}]) == i).^2)./abserr97(month([K_time{:}]) == i)));
    end
end

%% Separating summer and winter
H_hpss_summer = zeros(1, length(H_noy));
H_hpss_winter = zeros(1, length(H_noy));
K_hpss_summer = zeros(1, length(K_noy));
K_hpss_winter = zeros(1, length(K_noy));
H_sdep_summer = zeros(1, length(H_noy));
H_sdep_winter = zeros(1, length(H_noy));
K_sdep_summer = zeros(1, length(K_noy));
K_sdep_winter = zeros(1, length(K_noy));
H_LST_winter = cell(1, length(H_noy));
H_LST_winter(1,:) = {0};
H_LST_summer = cell(1, length(H_noy));
H_LST_summer(1,:) = {0};
K_LST_winter = cell(1, length(K_noy));
K_LST_winter(1,:) = {0};
K_LST_summer = cell(1, length(K_noy));
K_LST_summer(1,:) = {0};

intensity31_summer = zeros(1, length(H_noy));
intensity31_winter = zeros(1, length(H_noy));
intensity42_summer = zeros(1, length(H_noy));
intensity42_winter = zeros(1, length(H_noy));
intensity53_summer = zeros(1, length(H_noy));
intensity53_winter = zeros(1, length(H_noy));
intensity64_summer = zeros(1, length(H_noy));
intensity64_winter = zeros(1, length(H_noy));
intensity97_summer = zeros(1, length(K_noy));
intensity97_winter = zeros(1, length(K_noy));

j = 0;
k = 0;
for i = 1:length(H_noy)
    if H_noy(i) >= 79 && H_noy(i) <= 266
        j = j + 1;
        H_hpss_summer(j) = H_hpss(i);
        H_sdep_summer(j) = H_sdep(i);
        H_LST_summer{j} = H_LST{i};
        intensity31_summer(j) = intensity31(i);
        intensity42_summer(j) = intensity42(i);
        intensity53_summer(j) = intensity53(i);
        intensity64_summer(j) = intensity64(i);
    else
        k = k + 1;
        H_hpss_winter(k) = H_hpss(i);
        H_sdep_winter(k) = H_sdep(i);
        H_LST_winter{k} = H_LST{i};
        intensity31_winter(k) = intensity31(i);
        intensity42_winter(k) = intensity42(i);
        intensity53_winter(k) = intensity53(i);
        intensity64_winter(k) = intensity64(i);
    end
end
H_hpss_summer = H_hpss_summer(1:j);
H_sdep_summer = H_sdep_summer(1:j);
H_LST_summer = H_LST_summer(1:j);
intensity31_summer = intensity31_summer(1:j);
intensity42_summer = intensity42_summer(1:j);
intensity53_summer = intensity53_summer(1:j);
intensity64_summer = intensity64_summer(1:j);
H_hpss_winter = H_hpss_winter(1:k);
H_sdep_winter = H_sdep_winter(1:k);
H_LST_winter = H_LST_winter(1:k);
intensity31_winter = intensity31_winter(1:k);
intensity42_winter = intensity42_winter(1:k);
intensity53_winter = intensity53_winter(1:k);
intensity64_winter = intensity64_winter(1:k);

j = 0;
k = 0;
for i = 1:length(K_noy)
    if K_noy(i) >= 79 && K_noy(i) <= 266
        j = j + 1;
        K_hpss_summer(j) = K_hpss(i);
        K_sdep_summer(j) = K_sdep(i);
        K_LST_summer{j} = K_LST{i};
        intensity97_summer(j) = intensity97(i);
    else
        k = k + 1;
        K_hpss_winter(k) = K_hpss(i);
        K_sdep_winter(k) = K_sdep(i);
        K_LST_winter{k} = K_LST{i};
        intensity97_winter(k) = intensity97(i);
    end
end
K_hpss_summer = K_hpss_summer(1:j);
K_sdep_summer = K_sdep_summer(1:j);
K_LST_summer = K_LST_summer(1:j);
intensity97_summer = intensity97_summer(1:j);
K_hpss_winter = K_hpss_winter(1:k);
K_sdep_winter = K_sdep_winter(1:k);
K_LST_winter = K_LST_winter(1:k);
intensity97_winter = intensity97_winter(1:k);

%% Hours after sunset
%% Average of every hour
%% Combining seasons
H_hPostSunSet = hPostSunSet; % A quick fix for consistensy

intensity31_hpss_lowres = zeros(1,14);
intensity42_hpss_lowres = zeros(1,14);
intensity53_hpss_lowres = zeros(1,14);
intensity64_hpss_lowres = zeros(1,14);
intensity97_hpss_lowres = zeros(1,14);

intensity31_hpss_lowres_err = zeros(1,14);
intensity42_hpss_lowres_err = zeros(1,14);
intensity53_hpss_lowres_err = zeros(1,14);
intensity64_hpss_lowres_err = zeros(1,14);
intensity97_hpss_lowres_err = zeros(1,14);
%% Separating winter and summer init hpss
intensity31_hpss_lowres_summer = zeros(1,14);
intensity42_hpss_lowres_summer = zeros(1,14);
intensity53_hpss_lowres_summer = zeros(1,14);
intensity64_hpss_lowres_summer = zeros(1,14);
intensity97_hpss_lowres_summer = zeros(1,14);
intensity31_hpss_lowres_winter = zeros(1,14);
intensity42_hpss_lowres_winter = zeros(1,14);
intensity53_hpss_lowres_winter = zeros(1,14);
intensity64_hpss_lowres_winter = zeros(1,14);
intensity97_hpss_lowres_winter = zeros(1,14);

intensity31_hpss_lowres_summer_err = zeros(1,14);
intensity42_hpss_lowres_summer_err = zeros(1,14);
intensity53_hpss_lowres_summer_err = zeros(1,14);
intensity64_hpss_lowres_summer_err = zeros(1,14);
intensity97_hpss_lowres_summer_err = zeros(1,14);
intensity31_hpss_lowres_winter_err = zeros(1,14);
intensity42_hpss_lowres_winter_err = zeros(1,14);
intensity53_hpss_lowres_winter_err = zeros(1,14);
intensity64_hpss_lowres_winter_err = zeros(1,14);
intensity97_hpss_lowres_winter_err = zeros(1,14);

%% Common for-loop for data processing
for i = 0:13
    intensity31_hpss_lowres(i+1) = mean(intensity31(floor(H_hPostSunSet) == i));
    intensity42_hpss_lowres(i+1) = mean(intensity42(floor(H_hPostSunSet) == i));
    intensity53_hpss_lowres(i+1) = mean(intensity53(floor(H_hPostSunSet) == i));
    intensity64_hpss_lowres(i+1) = mean(intensity64(floor(H_hPostSunSet) == i));
    intensity97_hpss_lowres(i+1) = mean(intensity97(floor(K_hPostSunSet) == i));
    
    intensity31_hpss_lowres_err(i+1) = std(intensity31(floor(H_hPostSunSet) == i))...
        /length(intensity31(floor(H_hPostSunSet) == i));
    intensity42_hpss_lowres_err(i+1) = std(intensity42(floor(H_hPostSunSet) == i))...
        /length(intensity42(floor(H_hPostSunSet) == i));
    intensity53_hpss_lowres_err(i+1) = std(intensity53(floor(H_hPostSunSet) == i))...
        /length(intensity53(floor(H_hPostSunSet) == i));
    intensity64_hpss_lowres_err(i+1) = std(intensity64(floor(H_hPostSunSet) == i))...
        /length(intensity64(floor(H_hPostSunSet) == i));
    intensity97_hpss_lowres_err(i+1) = std(intensity97(floor(H_hPostSunSet) == i))...
        /length(intensity97(floor(K_hPostSunSet) == i));
end

%% Average of half hours, combining seasons
intensity31_hpss_highres = zeros(1,28);
intensity42_hpss_highres = zeros(1,28);
intensity53_hpss_highres = zeros(1,28);
intensity64_hpss_highres = zeros(1,28);
intensity97_hpss_highres = zeros(1,28);

intensity31_hpss_highres_err = zeros(1,28);
intensity42_hpss_highres_err = zeros(1,28);
intensity53_hpss_highres_err = zeros(1,28);
intensity64_hpss_highres_err = zeros(1,28);
intensity97_hpss_highres_err = zeros(1,28);

for i = 0:27
    intensity31_hpss_highres(i+1) = mean(intensity31(floor(H_hPostSunSet*2) == i));
    intensity42_hpss_highres(i+1) = mean(intensity42(floor(H_hPostSunSet*2) == i));
    intensity53_hpss_highres(i+1) = mean(intensity53(floor(H_hPostSunSet*2) == i));
    intensity64_hpss_highres(i+1) = mean(intensity64(floor(H_hPostSunSet*2) == i));
    intensity97_hpss_highres(i+1) = mean(intensity97(floor(H_hPostSunSet*2) == i));
    
    intensity31_hpss_highres_err(i+1) = std(intensity31(floor(H_hPostSunSet*2) == i))...
        /length(intensity31(floor(H_hPostSunSet*2) == i));
    intensity42_hpss_highres_err(i+1) = std(intensity42(floor(H_hPostSunSet*2) == i))...
        /length(intensity42(floor(H_hPostSunSet*2) == i));
    intensity53_hpss_highres_err(i+1) = std(intensity53(floor(H_hPostSunSet*2) == i))...
        /length(intensity53(floor(H_hPostSunSet*2) == i));
    intensity64_hpss_highres_err(i+1) = std(intensity64(floor(H_hPostSunSet*2) == i))...
        /length(intensity31(floor(H_hPostSunSet*2) == i));
    intensity97_hpss_highres_err(i+1) = std(intensity97(floor(H_hPostSunSet*2) == i))...
        /length(intensity97(floor(H_hPostSunSet*2) == i));
end

%% Eye ball plots
% Monthly averages of intensity vs hours after sunset/local time. Combine
% all years

eyeball97_moy_lowres = zeros(14,13); 
eyeball97_moy_highres = zeros(28,13);
eyeball97_noy_lowres = zeros(14, 366);
eyeball97_noy_highres = zeros(28,366);

%% Montly-hourly mean

for m = 1:12
    intensity97_m = intensity97(month([K_time{:}]) == m);
    K_hpss_m = K_hpss(month([K_time{:}]) == m);
    for h = 0:12  
        eyeball97_moy_lowres(h+1, m) = mean(intensity97_m(floor(K_hpss_m) == h));
    end
end

%% Create interpolation
interp_method = 2; % 2 or 4 (or 0,2
eyeball97_moy_lowres_interp_full = inpaint_nans(repmat(eyeball97_moy_lowres([1:13],[1:12]),3,3), interp_method);
eyeball97_moy_lowres_interp = eyeball97_moy_lowres_interp_full([14:27], [13:25]);

%% Monthly-half-hourly mean
for m = 1:12
    intensity97_m = intensity97(month([K_time{:}]) == m);
    K_hpss_m = K_hpss(month([K_time{:}]) == m);
    for h = 0:24  
        eyeball97_moy_highres(h+1, m) = mean(intensity97_m(floor(K_hpss_m*2) == h));
    end
end

%% Create interpolation
eyeball97_moy_highres_interp_full = inpaint_nans(repmat(eyeball97_moy_highres([1:25],[1:12]),3,3), interp_method); % Use method 2 or 4
eyeball97_moy_highres_interp = eyeball97_moy_highres_interp_full([26:51], [13:25]);

%% Night of year-hourly mean
for n = 1:365
    intensity97_n = intensity97(K_noy == n);
    K_hpss_n = K_hpss(K_noy == n);
    for h = 0:12  
        eyeball97_noy_lowres(h+1, n) = mean(intensity97_n(floor(K_hpss_n) == h));
    end
end

%% Create interpolation
eyeball97_noy_lowres_interp_full = inpaint_nans(repmat(eyeball97_noy_lowres([1:13],[1:365]),3,3), interp_method);
eyeball97_noy_lowres_interp = eyeball97_noy_lowres_interp_full([14:27], [366:731]);

%% Night of year-half-hourly mean

for n = 1:365
    intensity97_n = intensity97(K_noy == n);
    K_hpss_n = K_hpss(K_noy == n);
    for h = 0:24  
        eyeball97_noy_highres(h+1, n) = mean(intensity97_n(floor(K_hpss_n*2) == h));
    end
end

%% Create interpolation
eyeball97_noy_highres_interp_full = inpaint_nans(repmat(eyeball97_noy_highres([1:25],[1:365]),3,3), interp_method);
eyeball97_noy_highres_interp = eyeball97_noy_highres_interp_full([26:51], [366:731]);

%% Solar depression angles
%% Low resolution
%% Variable init
N = 18;
midnight = 6.0; % Hours after sunset
midnight_summer = 5.5;
midnight_winter = 6.5;

intensity31_setting = intensity31(hPostSunSet <= midnight);
intensity31_rising = intensity31(hPostSunSet >= midnight); % some overlap
intensity42_setting = intensity42(hPostSunSet <= midnight);
intensity42_rising = intensity42(hPostSunSet >= midnight); % some overlap
intensity53_setting = intensity53(hPostSunSet <= midnight);
intensity53_rising = intensity53(hPostSunSet >= midnight); % some overlap
intensity64_setting = intensity64(hPostSunSet <= midnight);
intensity64_rising = intensity64(hPostSunSet >= midnight); % some overlap

intensity31_summer_setting = intensity31_summer(H_hpss_summer <= midnight_summer);
intensity31_summer_rising = intensity31_summer(H_hpss_summer >= midnight_summer);
intensity31_winter_setting = intensity31_winter(H_hpss_winter <= midnight_winter);
intensity31_winter_rising = intensity31_winter(H_hpss_winter >= midnight_winter);
intensity42_summer_setting = intensity42_summer(H_hpss_summer <= midnight_summer);
intensity42_summer_rising = intensity42_summer(H_hpss_summer >= midnight_summer);
intensity42_winter_setting = intensity42_winter(H_hpss_winter <= midnight_winter);
intensity42_winter_rising = intensity42_winter(H_hpss_winter >= midnight_winter);
intensity53_summer_setting = intensity53_summer(H_hpss_summer <= midnight_summer);
intensity53_summer_rising = intensity53_summer(H_hpss_summer >= midnight_summer);
intensity53_winter_setting = intensity53_winter(H_hpss_winter <= midnight_winter);
intensity53_winter_rising = intensity53_winter(H_hpss_winter >= midnight_winter);
intensity64_summer_setting = intensity64_summer(H_hpss_summer <= midnight_summer);
intensity64_summer_rising = intensity64_summer(H_hpss_summer >= midnight_summer);
intensity64_winter_setting = intensity64_winter(H_hpss_winter <= midnight_winter);
intensity64_winter_rising = intensity64_winter(H_hpss_winter >= midnight_winter);

H_sdep_setting = H_sdep(hPostSunSet <= midnight);
H_sdep_rising = H_sdep(hPostSunSet >= midnight);

H_sdep_summer_setting = H_sdep_summer(H_hpss_summer <= midnight_summer);
H_sdep_summer_rising = H_sdep_summer(H_hpss_summer >= midnight_summer);
H_sdep_winter_setting = H_sdep_winter(H_hpss_winter <= midnight_winter);
H_sdep_winter_rising = H_sdep_winter(H_hpss_winter >= midnight_winter);

intensity97_setting = intensity97(K_hPostSunSet <= midnight);
intensity97_rising = intensity97(K_hPostSunSet >= midnight); % some overlap

intensity97_summer_setting = intensity97_summer(K_hpss_summer <= midnight_summer);
intensity97_summer_rising = intensity97_summer(K_hpss_summer >= midnight_summer);
intensity97_winter_setting = intensity97_winter(K_hpss_winter <= midnight_winter);
intensity97_winter_rising = intensity97_winter(K_hpss_winter >= midnight_winter);

K_sdep_setting = K_sdep(K_hPostSunSet <= midnight);
K_sdep_rising = K_sdep(K_hPostSunSet >= midnight);
K_sdep_summer_setting = K_sdep_summer(K_hpss_summer <= midnight_summer);
K_sdep_summer_rising = K_sdep_summer(K_hpss_summer >= midnight_summer);
K_sdep_winter_setting = K_sdep_winter(K_hpss_winter <= midnight_winter);
K_sdep_winter_rising = K_sdep_winter(K_hpss_winter >= midnight_winter);

K_time_setting = K_time(K_hPostSunSet <= midnight);
K_time_rising = K_time(K_hPostSunSet >= midnight);
H_time_setting = time(H_hPostSunSet <= midnight);
H_time_rising = time(H_hPostSunSet >= midnight);

intensity31_sdep_setting_lowres = zeros(1,N);
intensity31_sdep_rising_lowres = zeros(1,N);
intensity42_sdep_setting_lowres = zeros(1,N);
intensity42_sdep_rising_lowres = zeros(1,N);
intensity53_sdep_setting_lowres = zeros(1,N);
intensity53_sdep_rising_lowres = zeros(1,N);
intensity64_sdep_setting_lowres = zeros(1,N);
intensity64_sdep_rising_lowres = zeros(1,N);
intensity97_sdep_setting_lowres = zeros(1,N);
intensity97_sdep_rising_lowres = zeros(1,N);

intensity31_sdep_summer_setting_lowres = zeros(1,N);
intensity31_sdep_summer_rising_lowres = zeros(1,N);
intensity42_sdep_summer_setting_lowres = zeros(1,N);
intensity42_sdep_summer_rising_lowres = zeros(1,N);
intensity53_sdep_summer_setting_lowres = zeros(1,N);
intensity53_sdep_summer_rising_lowres = zeros(1,N);
intensity64_sdep_summer_setting_lowres = zeros(1,N);
intensity64_sdep_summer_rising_lowres = zeros(1,N);
intensity97_sdep_summer_setting_lowres = zeros(1,N);
intensity97_sdep_summer_rising_lowres = zeros(1,N);
intensity31_sdep_winter_setting_lowres = zeros(1,N);
intensity31_sdep_winter_rising_lowres = zeros(1,N);
intensity42_sdep_winter_setting_lowres = zeros(1,N);
intensity42_sdep_winter_rising_lowres = zeros(1,N);
intensity53_sdep_winter_setting_lowres = zeros(1,N);
intensity53_sdep_winter_rising_lowres = zeros(1,N);
intensity64_sdep_winter_setting_lowres = zeros(1,N);
intensity64_sdep_winter_rising_lowres = zeros(1,N);
intensity97_sdep_winter_setting_lowres = zeros(1,N);
intensity97_sdep_winter_rising_lowres = zeros(1,N);

intensity31_sdep_setting_lowres_err = zeros(1,N);
intensity31_sdep_rising_lowres_err = zeros(1,N);
intensity42_sdep_setting_lowres_err = zeros(1,N);
intensity42_sdep_rising_lowres_err = zeros(1,N);
intensity53_sdep_setting_lowres_err = zeros(1,N);
intensity53_sdep_rising_lowres_err = zeros(1,N);
intensity64_sdep_setting_lowres_err = zeros(1,N);
intensity64_sdep_rising_lowres_err = zeros(1,N);
intensity97_sdep_setting_lowres_err = zeros(1,N);
intensity97_sdep_rising_lowres_err = zeros(1,N);
intensity31_sdep_summer_setting_lowres_err = zeros(1,N);
intensity31_sdep_summer_rising_lowres_err = zeros(1,N);
intensity42_sdep_summer_setting_lowres_err = zeros(1,N);
intensity42_sdep_summer_rising_lowres_err = zeros(1,N);
intensity53_sdep_summer_setting_lowres_err = zeros(1,N);
intensity53_sdep_summer_rising_lowres_err = zeros(1,N);
intensity64_sdep_summer_setting_lowres_err = zeros(1,N);
intensity64_sdep_summer_rising_lowres_err = zeros(1,N);
intensity97_sdep_summer_setting_lowres_err = zeros(1,N);
intensity97_sdep_summer_rising_lowres_err = zeros(1,N);
intensity31_sdep_winter_setting_lowres_err = zeros(1,N);
intensity31_sdep_winter_rising_lowres_err = zeros(1,N);
intensity42_sdep_winter_setting_lowres_err = zeros(1,N);
intensity42_sdep_winter_rising_lowres_err = zeros(1,N);
intensity53_sdep_winter_setting_lowres_err = zeros(1,N);
intensity53_sdep_winter_rising_lowres_err = zeros(1,N);
intensity64_sdep_winter_setting_lowres_err = zeros(1,N);
intensity64_sdep_winter_rising_lowres_err = zeros(1,N);
intensity97_sdep_winter_setting_lowres_err = zeros(1,N);
intensity97_sdep_winter_rising_lowres_err = zeros(1,N);

%% For loop data processing
for i = 1:18
    %% Combining seasons
    %% Setting
    temp_int97 = zeros(1, length(K_sdep_setting));
    temp_int31 = zeros(1, length(H_sdep_setting));
    temp_int42 = zeros(1, length(H_sdep_setting));
    temp_int53 = zeros(1, length(H_sdep_setting));
    temp_int64 = zeros(1, length(H_sdep_setting));
    k = 0;
    h = 0;
    for j = 1:length(K_sdep_setting)
        if (K_sdep_setting(j) > 0 + 5*(i-1)) && (K_sdep_setting(j) < 5 + 5*(i-1))
            k = k + 1;
            temp_int97(k) = intensity97_setting(j);
        end
    end
    for j = 1:length(H_sdep_setting)
        if (H_sdep_setting(j) > 0 + 5*(i-1)) && (H_sdep_setting(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_setting(j);
            temp_int42(h) = intensity42_setting(j);
            temp_int53(h) = intensity53_setting(j);
            temp_int64(h) = intensity64_setting(j);
        end
    end

    temp_int97 = temp_int97(1:k); % Remove zero-tail
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    intensity31_sdep_setting_lowres(i) = mean(temp_int31);
    intensity42_sdep_setting_lowres(i) = mean(temp_int42);
    intensity53_sdep_setting_lowres(i) = mean(temp_int53);
    intensity64_sdep_setting_lowres(i) = mean(temp_int64);
    intensity97_sdep_setting_lowres(i) = mean(temp_int97);
    intensity31_sdep_setting_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_setting_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_setting_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_setting_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_setting_lowres_err(i) = std(temp_int97)/length(temp_int97);
    clearvars temp_int*; 
    %% Rising

    temp_int97 = zeros(1, length(K_sdep_rising));
    temp_int31 = zeros(1, length(H_sdep_rising));
    temp_int42 = zeros(1, length(H_sdep_rising));
    temp_int53 = zeros(1, length(H_sdep_rising));
    temp_int64 = zeros(1, length(H_sdep_rising));
    k = 0;
    h = 0;
    for j = 1:length(K_sdep_rising)
        if (K_sdep_rising(j) > 0 + 5*(i-1)) && (K_sdep_rising(j) < 5 + 5*(i-1))
            k = k + 1;
            temp_int97(k) = intensity97_rising(j);
        end
    end
    for j = 1:length(H_sdep_rising)
        if (H_sdep_rising(j) > 0 + 5*(i-1)) && (H_sdep_rising(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_rising(j);
            temp_int42(h) = intensity42_rising(j);
            temp_int53(h) = intensity53_rising(j);
            temp_int64(h) = intensity64_rising(j);
        end
    end

    temp_int97 = temp_int97(1:k); % Remove zero-tail
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    intensity31_sdep_rising_lowres(i) = mean(temp_int31);
    intensity42_sdep_rising_lowres(i) = mean(temp_int42);
    intensity53_sdep_rising_lowres(i) = mean(temp_int53);
    intensity64_sdep_rising_lowres(i) = mean(temp_int64);
    intensity97_sdep_rising_lowres(i) = mean(temp_int97);
    intensity31_sdep_rising_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_rising_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_rising_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_rising_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_rising_lowres_err(i) = std(temp_int97)/length(temp_int97);
    clearvars temp_int*;
    %% Winter - summer
    %% Winter | Setting
    temp_int31 = zeros(1, length(H_sdep_winter_setting));
    temp_int42 = zeros(1, length(H_sdep_winter_setting));
    temp_int53 = zeros(1, length(H_sdep_winter_setting));
    temp_int64 = zeros(1, length(H_sdep_winter_setting));
    temp_int97 = zeros(1, length(K_sdep_winter_setting));

    k = 0;
    h = 0;
    for j = 1:length(H_sdep_winter_setting)
        if (H_sdep_winter_setting(j) > 0 + 5*(i-1)) && (H_sdep_winter_setting(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_winter_setting(j);
            temp_int42(h) = intensity42_winter_setting(j);
            temp_int53(h) = intensity53_winter_setting(j);
            temp_int64(h) = intensity64_winter_setting(j);
        end
    end
    for j = 1:length(K_sdep_winter_setting)
        if (K_sdep_winter_setting(j) > 0 + 5*(i-1)) && (K_sdep_winter_setting(j) < 5 + 5*(i-1))
            %%
            k = k + 1;
            temp_int97(k) = intensity97_winter_setting(j);
        end
    end
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    temp_int97 = temp_int97(1:k);
    intensity31_sdep_winter_setting_lowres(i) = mean(temp_int31);
    intensity42_sdep_winter_setting_lowres(i) = mean(temp_int42);
    intensity53_sdep_winter_setting_lowres(i) = mean(temp_int53);
    intensity64_sdep_winter_setting_lowres(i) = mean(temp_int64);
    intensity97_sdep_winter_setting_lowres(i) = mean(temp_int97);
    intensity31_sdep_winter_setting_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_winter_setting_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_winter_setting_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_winter_setting_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_winter_setting_lowres_err(i) = std(temp_int97)/length(temp_int97);
    
    clearvars h k temp_int*
    %% Winter | Rising
    temp_int31 = zeros(1, length(H_sdep_winter_rising));
    temp_int42 = zeros(1, length(H_sdep_winter_rising));
    temp_int53 = zeros(1, length(H_sdep_winter_rising));
    temp_int64 = zeros(1, length(H_sdep_winter_rising));
    temp_int97 = zeros(1, length(K_sdep_winter_rising));

    k = 0;
    h = 0;
    for j = 1:length(H_sdep_winter_rising)
        if (H_sdep_winter_rising(j) > 0 + 5*(i-1)) && (H_sdep_winter_rising(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_winter_rising(j);
            temp_int42(h) = intensity42_winter_rising(j);
            temp_int53(h) = intensity53_winter_rising(j);
            temp_int64(h) = intensity64_winter_rising(j);
        end
    end
    for j = 1:length(K_sdep_winter_rising)
        if (K_sdep_winter_rising(j) > 0 + 5*(i-1)) && (K_sdep_winter_rising(j) < 5 + 5*(i-1))
            %%
            k = k + 1;
            temp_int97(k) = intensity97_winter_rising(j);
        end
    end
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    temp_int97 = temp_int97(1:k);
    intensity31_sdep_winter_rising_lowres(i) = mean(temp_int31);
    intensity42_sdep_winter_rising_lowres(i) = mean(temp_int42);
    intensity53_sdep_winter_rising_lowres(i) = mean(temp_int53);
    intensity64_sdep_winter_rising_lowres(i) = mean(temp_int64);
    intensity97_sdep_winter_rising_lowres(i) = mean(temp_int97);
    intensity31_sdep_winter_rising_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_winter_rising_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_winter_rising_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_winter_rising_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_winter_rising_lowres_err(i) = std(temp_int97)/length(temp_int97);
    
    clearvars h k temp_int*
    %% Summer | setting
    temp_int31 = zeros(1, length(H_sdep_summer_setting));
    temp_int42 = zeros(1, length(H_sdep_summer_setting));
    temp_int53 = zeros(1, length(H_sdep_summer_setting));
    temp_int64 = zeros(1, length(H_sdep_summer_setting));
    temp_int97 = zeros(1, length(K_sdep_summer_setting));

    k = 0;
    h = 0;
    for j = 1:length(H_sdep_summer_setting)
        if (H_sdep_summer_setting(j) > 0 + 5*(i-1)) && (H_sdep_summer_setting(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_summer_setting(j);
            temp_int42(h) = intensity42_summer_setting(j);
            temp_int53(h) = intensity53_summer_setting(j);
            temp_int64(h) = intensity64_summer_setting(j);
        end
    end
    for j = 1:length(K_sdep_summer_setting)
        if (K_sdep_summer_setting(j) > 0 + 5*(i-1)) && (K_sdep_summer_setting(j) < 5 + 5*(i-1))
            %%
            k = k + 1;
            temp_int97(k) = intensity97_summer_setting(j);
        end
    end
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    temp_int97 = temp_int97(1:k);
    intensity31_sdep_summer_setting_lowres(i) = mean(temp_int31);
    intensity42_sdep_summer_setting_lowres(i) = mean(temp_int42);
    intensity53_sdep_summer_setting_lowres(i) = mean(temp_int53);
    intensity64_sdep_summer_setting_lowres(i) = mean(temp_int64);
    intensity97_sdep_summer_setting_lowres(i) = mean(temp_int97);
    intensity31_sdep_summer_setting_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_summer_setting_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_summer_setting_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_summer_setting_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_summer_setting_lowres_err(i) = std(temp_int97)/length(temp_int97);
    
    clearvars h k temp_int*
    %% Summer | rising
    temp_int31 = zeros(1, length(H_sdep_summer_rising));
    temp_int42 = zeros(1, length(H_sdep_summer_rising));
    temp_int53 = zeros(1, length(H_sdep_summer_rising));
    temp_int64 = zeros(1, length(H_sdep_summer_rising));
    temp_int97 = zeros(1, length(K_sdep_summer_rising));

    k = 0;
    h = 0;
    for j = 1:length(H_sdep_summer_rising)
        if (H_sdep_summer_rising(j) > 0 + 5*(i-1)) && (H_sdep_summer_rising(j) < 5 + 5*(i-1))
            %%
            h = h + 1;
            temp_int31(h) = intensity31_summer_rising(j);
            temp_int42(h) = intensity42_summer_rising(j);
            temp_int53(h) = intensity53_summer_rising(j);
            temp_int64(h) = intensity64_summer_rising(j);
        end
    end
    for j = 1:length(K_sdep_summer_rising)
        if (K_sdep_summer_rising(j) > 0 + 5*(i-1)) && (K_sdep_summer_rising(j) < 5 + 5*(i-1))
            %%
            k = k + 1;
            temp_int97(k) = intensity97_summer_rising(j);
        end
    end
    temp_int31 = temp_int31(1:h);
    temp_int42 = temp_int42(1:h);
    temp_int53 = temp_int53(1:h);
    temp_int64 = temp_int64(1:h);
    temp_int97 = temp_int97(1:k);
    intensity31_sdep_summer_rising_lowres(i) = mean(temp_int31);
    intensity42_sdep_summer_rising_lowres(i) = mean(temp_int42);
    intensity53_sdep_summer_rising_lowres(i) = mean(temp_int53);
    intensity64_sdep_summer_rising_lowres(i) = mean(temp_int64);
    intensity97_sdep_summer_rising_lowres(i) = mean(temp_int97);
    intensity31_sdep_summer_rising_lowres_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_sdep_summer_rising_lowres_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_sdep_summer_rising_lowres_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_sdep_summer_rising_lowres_err(i) = std(temp_int64)/length(temp_int64);
    intensity97_sdep_summer_rising_lowres_err(i) = std(temp_int97)/length(temp_int97);
    
    clearvars h k temp_int*
end

%% Solar depression angle eyeball plot
eyeball_sdep_97_setting = zeros(18, 13);
eyeball_sdep_97_rising = zeros(18, 13);

for m = 1:12
    for i = 1:18
        K_sdep_setting_m = K_sdep_setting(month([K_time_setting{:}]) == m); % Temporary variables for each month
        K_sdep_rising_m = K_sdep_rising(month([K_time_rising{:}]) == m);
        intensity97_setting_m = intensity97_setting(month([K_time_setting{:}]) == m);
        intensity97_rising_m = intensity97_rising(month([K_time_rising{:}]) == m);
        
        %% Setting
        temp_int97 = zeros(1, length(K_sdep_setting_m));
        k = 0;
        for j = 1:length(K_sdep_setting_m)
            if K_sdep_setting_m(j) >= 0 + 5*(i-1) && K_sdep_setting_m(j) <= 5 + 5*(i-1)
                k = k + 1;
                temp_int97(k) = intensity97_setting_m(j);
            end
        end
        temp_int97 = temp_int97(1:k); % remove zero-tail
        
        eyeball_sdep_97_setting(i, m) = mean(temp_int97);
        clearvars k temp_int97;
        
        %% Rising
        temp_int97 = zeros(1, length(K_sdep_rising_m));
        k = 0;
        for j = 1:length(K_sdep_rising_m)
            if K_sdep_rising_m(j) >= 0 + 5*(i-1) && K_sdep_rising_m(j) <= 5 + 5*(i-1)
                k = k + 1;
                temp_int97(k) = intensity97_rising_m(j);
            end
        end
        temp_int97 = temp_int97(1:k); % remove zero-tail
        
        eyeball_sdep_97_rising(i, m) = mean(temp_int97);
        clearvars k temp_int97;
    end
end

eyeball_sdep_97 = [eyeball_sdep_97_setting; flipud(eyeball_sdep_97_rising)];

%% High Resolution
N = 36; % Higher resolution
intensity31_sdep = zeros(1,N);
intensity31_bm = intensity31(hour([time{:}]) > 15);
intensity31_am = intensity31(hour([time{:}]) < 10);
H_sdep_bm = H_sdep(hour([time{:}]) > 15);
H_sdep_am = H_sdep(hour([time{:}]) < 10);

%% Setting - rising, Combining seasons

midnight = 6.0; % Hours after sunset
intensity31_setting = intensity31(hPostSunSet <= midnight);
intensity31_rising = intensity31(hPostSunSet >= midnight); % some overlap
intensity42_setting = intensity42(hPostSunSet <= midnight);
intensity42_rising = intensity42(hPostSunSet >= midnight); % some overlap
intensity53_setting = intensity53(hPostSunSet <= midnight);
intensity53_rising = intensity53(hPostSunSet >= midnight); % some overlap
intensity64_setting = intensity64(hPostSunSet <= midnight);
intensity64_rising = intensity64(hPostSunSet >= midnight); % some overlap

H_sdep_setting = H_sdep(hPostSunSet <= midnight);
H_sdep_rising = H_sdep(hPostSunSet >= midnight);

intensity97_setting = intensity97(K_hPostSunSet <= midnight);
intensity97_rising = intensity97(K_hPostSunSet >= midnight); % some overlap

K_time_setting = K_time(K_hPostSunSet <= midnight);
K_time_rising = K_time(K_hPostSunSet >= midnight);
H_time_setting = time(H_hPostSunSet <= midnight);
H_time_rising = time(H_hPostSunSet >= midnight);


K_sdep_setting = K_sdep(K_hPostSunSet <= midnight);
K_sdep_rising = K_sdep(K_hPostSunSet >= midnight);

intensity31_sdep_setting = zeros(1,N);
intensity31_sdep_rising = zeros(1,N);
intensity42_sdep_setting = zeros(1,N);
intensity42_sdep_rising = zeros(1,N);
intensity53_sdep_setting = zeros(1,N);
intensity53_sdep_rising = zeros(1,N);
intensity64_sdep_setting = zeros(1,N);
intensity64_sdep_rising = zeros(1,N);
intensity97_sdep_setting = zeros(1,N);
intensity97_sdep_rising = zeros(1,N);

intensity31_sdep_setting_err = zeros(1,N);
intensity31_sdep_rising_err = zeros(1,N);
intensity42_sdep_setting_err = zeros(1,N);
intensity42_sdep_rising_err = zeros(1,N);
intensity53_sdep_setting_err = zeros(1,N);
intensity53_sdep_rising_err = zeros(1,N);
intensity64_sdep_setting_err = zeros(1,N);
intensity64_sdep_rising_err = zeros(1,N);
intensity97_sdep_setting_err = zeros(1,N);
intensity97_sdep_rising_err = zeros(1,N);

%% Setting - rising | Summer - winter

intensity31_summer_setting = intensity31_summer(H_hpss_summer <= 5.5);
intensity31_summer_rising = intensity31_summer(H_hpss_summer >= 5.5);
intensity42_summer_setting = intensity42_summer(H_hpss_summer <= 5.5);
intensity42_summer_rising = intensity42_summer(H_hpss_summer >= 5.5);
intensity53_summer_setting = intensity53_summer(H_hpss_summer <= 5.5);
intensity53_summer_rising = intensity53_summer(H_hpss_summer >= 5.5);
intensity64_summer_setting = intensity53_summer(H_hpss_summer <= 5.5);
intensity64_summer_rising = intensity64_summer(H_hpss_summer >= 5.5);

intensity31_winter_setting = intensity31_winter(H_hpss_winter <= 6.5);
intensity31_winter_rising = intensity31_winter(H_hpss_winter >= 6.5);
intensity42_winter_setting = intensity42_winter(H_hpss_winter <= 6.5);
intensity42_winter_rising = intensity42_winter(H_hpss_winter >= 6.5);
intensity53_winter_setting = intensity53_winter(H_hpss_winter <= 6.5);
intensity53_winter_rising = intensity53_winter(H_hpss_winter >= 6.5);
intensity64_winter_setting = intensity53_winter(H_hpss_winter <= 6.5);
intensity64_winter_rising = intensity64_winter(H_hpss_winter >= 6.5);

H_sdep_summer_setting = H_sdep_summer(H_hpss_summer <= 5.5);
H_sdep_summer_rising = H_sdep_summer(H_hpss_summer >= 5.5);
H_sdep_winter_setting = H_sdep_winter(H_hpss_winter <= 6.5);
H_sdep_winter_rising = H_sdep_winter(H_hpss_winter >= 6.5);

intensity97_summer_setting = intensity97_summer(K_hpss_summer <= 5.5);
intensity97_summer_rising = intensity97_summer(K_hpss_summer >= 5.5);
intensity97_winter_setting = intensity97_winter(K_hpss_winter <= 6.5);
intensity97_winter_rising = intensity97_winter(K_hpss_winter >= 6.5);

K_sdep_summer_setting = K_sdep_summer(K_hpss_summer <= 5.5);
K_sdep_summer_rising = K_sdep_summer(K_hpss_summer >= 5.5);
K_sdep_winter_setting = K_sdep_winter(K_hpss_winter <= 6.5);
K_sdep_winter_rising = K_sdep_winter(K_hpss_winter >= 6.5);

intensity31_sdep_summer_setting = zeros(1,N);
intensity31_sdep_summer_rising = zeros(1,N);
intensity31_sdep_winter_setting = zeros(1,N);
intensity31_sdep_winter_rising = zeros(1,N);
intensity31_sdep_summer_setting_err = zeros(1,N);
intensity31_sdep_summer_rising_err = zeros(1,N);
intensity31_sdep_winter_setting_err = zeros(1,N);
intensity31_sdep_winter_rising_err = zeros(1,N);

intensity42_sdep_summer_setting = zeros(1,N);
intensity42_sdep_summer_rising = zeros(1,N);
intensity42_sdep_winter_setting = zeros(1,N);
intensity42_sdep_winter_rising = zeros(1,N);
intensity42_sdep_summer_setting_err = zeros(1,N);
intensity42_sdep_summer_rising_err = zeros(1,N);
intensity42_sdep_winter_setting_err = zeros(1,N);
intensity42_sdep_winter_rising_err = zeros(1,N);

intensity53_sdep_summer_setting = zeros(1,N);
intensity53_sdep_summer_rising = zeros(1,N);
intensity53_sdep_winter_setting = zeros(1,N);
intensity53_sdep_winter_rising = zeros(1,N);
intensity53_sdep_summer_setting_err = zeros(1,N);
intensity53_sdep_summer_rising_err = zeros(1,N);
intensity53_sdep_winter_setting_err = zeros(1,N);
intensity53_sdep_winter_rising_err = zeros(1,N);

intensity64_sdep_summer_setting = zeros(1,N);
intensity64_sdep_summer_rising = zeros(1,N);
intensity64_sdep_winter_setting = zeros(1,N);
intensity64_sdep_winter_rising = zeros(1,N);
intensity64_sdep_summer_setting_err = zeros(1,N);
intensity64_sdep_summer_rising_err = zeros(1,N);
intensity64_sdep_winter_setting_err = zeros(1,N);
intensity64_sdep_winter_rising_err = zeros(1,N);

intensity97_sdep_summer_setting = zeros(1,N);
intensity97_sdep_summer_rising = zeros(1,N);
intensity97_sdep_winter_setting = zeros(1,N);
intensity97_sdep_winter_rising = zeros(1,N);
intensity97_sdep_summer_setting_err = zeros(1,N);
intensity97_sdep_summer_rising_err = zeros(1,N);
intensity97_sdep_winter_setting_err = zeros(1,N);
intensity97_sdep_winter_rising_err = zeros(1,N);

intensity42_sdep = zeros(1,N);
intensity53_sdep = zeros(1,N);
intensity64_sdep = zeros(1,N);
intensity97_sdep = zeros(1,N);
intensity31_sdep_err = zeros(1,N);
intensity42_sdep_err = zeros(1,N);
intensity53_sdep_err = zeros(1,N);
intensity64_sdep_err = zeros(1,N);
intensity97_sdep_err = zeros(1,N);

for i = 1:N
    intensity97_sdep(i) = mean(intensity97(K_sdep(K_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    intensity97_sdep_err(i) = std(intensity97(K_sdep(K_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)))...
        /length(intensity97(K_sdep(K_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));

    intensity31_sdep(i) = mean(intensity31(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    intensity31_sdep_err(i) = std(intensity31(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)))...
        /length(intensity31(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    
    intensity42_sdep(i) = mean(intensity42(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    intensity42_sdep_err(i) = std(intensity42(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)))...
        /length(intensity97(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    
    intensity53_sdep(i) = mean(intensity53(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    intensity53_sdep_err(i) = std(intensity53(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)))...
        /length(intensity53(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    
    intensity64_sdep(i) = mean(intensity64(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    intensity64_sdep_err(i) = std(intensity64(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)))...
        /length(intensity64(H_sdep(H_sdep > 0 + 2.5*(i-1)) < 2.5 + 2.5*(i-1)));
    
    %% Combining seasons, setting-rising
    intensity31_sdep_setting(i) = mean(intensity31_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_rising(i) = mean(intensity31_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_setting(i) = mean(intensity42_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_rising(i) = mean(intensity42_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_setting(i) = mean(intensity53_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_rising(i) = mean(intensity53_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_setting(i) = mean(intensity64_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_rising(i) = mean(intensity64_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_setting(i) = mean(intensity97_setting(K_sdep_setting(K_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_rising(i) = mean(intensity97_rising(K_sdep_rising(K_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity31_sdep_setting_err(i) = std(intensity31_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_rising_err(i) = std(intensity31_rising(H_sdep_rising(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_setting_err(i) = std(intensity42_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_rising_err(i) = std(intensity42_rising(H_sdep_rising(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_setting_err(i) = std(intensity53_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_rising_err(i) = std(intensity53_rising(H_sdep_rising(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_setting_err(i) = std(intensity64_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_setting(H_sdep_setting(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_rising_err(i) = std(intensity64_rising(H_sdep_rising(H_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_rising(H_sdep_rising(H_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_setting_err(i) = std(intensity97_setting(K_sdep_setting(K_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_setting(K_sdep_setting(K_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_rising_err(i) = std(intensity97_rising(K_sdep_rising(K_sdep_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_rising(K_sdep_rising(K_sdep_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    %% Winter-summer setting-rising
    %% 31
    intensity31_sdep_winter_setting(i) = ...
        mean(intensity31_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_winter_setting_err(i) = ...
        std(intensity31_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity31_sdep_winter_rising(i) = ...
        mean(intensity31_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_winter_rising_err(i) = ...
        std(intensity31_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity31_sdep_summer_setting(i) = ...
        mean(intensity31_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_summer_setting_err(i) = ...
        std(intensity31_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity31_sdep_summer_rising(i) = ...
        mean(intensity31_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity31_sdep_summer_rising_err(i) = ...
        std(intensity31_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity31_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    %% 42
    intensity42_sdep_winter_setting(i) = ...
        mean(intensity42_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_winter_setting_err(i) = ...
        std(intensity42_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity42_sdep_winter_rising(i) = ...
        mean(intensity42_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_winter_rising_err(i) = ...
        std(intensity42_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity42_sdep_summer_setting(i) = ...
        mean(intensity42_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_summer_setting_err(i) = ...
        std(intensity42_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity42_sdep_summer_rising(i) = ...
        mean(intensity42_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity42_sdep_summer_rising_err(i) = ...
        std(intensity42_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity42_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1))); 
    %% 53
    intensity53_sdep_winter_setting(i) = ...
        mean(intensity53_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_winter_setting_err(i) = ...
        std(intensity53_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity53_sdep_winter_rising(i) = ...
        mean(intensity53_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_winter_rising_err(i) = ...
        std(intensity53_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity53_sdep_summer_setting(i) = ...
        mean(intensity53_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_summer_setting_err(i) = ...
        std(intensity53_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity53_sdep_summer_rising(i) = ...
        mean(intensity53_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity53_sdep_summer_rising_err(i) = ...
        std(intensity53_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity53_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1))); 
    %% 64
    intensity64_sdep_winter_setting(i) = ...
        mean(intensity64_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_winter_setting_err(i) = ...
        std(intensity64_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_winter_setting(H_sdep_winter_setting(H_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity64_sdep_winter_rising(i) = ...
        mean(intensity64_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_winter_rising_err(i) = ...
        std(intensity64_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_winter_rising(H_sdep_winter_rising(H_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity64_sdep_summer_setting(i) = ...
        mean(intensity64_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_summer_setting_err(i) = ...
        std(intensity64_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_summer_setting(H_sdep_summer_setting(H_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity64_sdep_summer_rising(i) = ...
        mean(intensity64_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity64_sdep_summer_rising_err(i) = ...
        std(intensity64_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity64_summer_rising(H_sdep_summer_rising(H_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1))); 
    %% 97 K filter
    intensity97_sdep_winter_setting(i) = ...
        mean(intensity97_winter_setting(K_sdep_winter_setting(K_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_winter_setting_err(i) = ...
        std(intensity97_winter_setting(K_sdep_winter_setting(K_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_winter_setting(K_sdep_winter_setting(K_sdep_winter_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity97_sdep_winter_rising(i) = ...
        mean(intensity97_winter_rising(K_sdep_winter_rising(K_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_winter_rising_err(i) = ...
        std(intensity97_winter_rising(K_sdep_winter_rising(K_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_winter_rising(K_sdep_winter_rising(K_sdep_winter_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity97_sdep_summer_setting(i) = ...
        mean(intensity97_summer_setting(K_sdep_summer_setting(K_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_summer_setting_err(i) = ...
        std(intensity97_summer_setting(K_sdep_summer_setting(K_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_summer_setting(K_sdep_summer_setting(K_sdep_summer_setting > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    
    intensity97_sdep_summer_rising(i) = ...
        mean(intensity97_summer_rising(K_sdep_summer_rising(K_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
    intensity97_sdep_summer_rising_err(i) = ...
        std(intensity97_summer_rising(K_sdep_summer_rising(K_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)))...
        /length(intensity97_summer_rising(K_sdep_summer_rising(K_sdep_summer_rising > 0 + 2.5*(i-1)) <= 2.5 + 2.5*(i-1)));
end

%% Local solar time
%% lowres
% We have LSTs from 18 to 23 and 0 to 5
H_lst_hours = zeros(1, length(H_LST));
K_lst_hours = zeros(1, length(K_LST));
lst_hours_lowres = [17 18 19 20 21 22 23 0 1 2 3 4 5 6];
intensity31_lst_lowres = zeros(1, 14);
intensity42_lst_lowres = zeros(1, 14);
intensity53_lst_lowres = zeros(1, 14);
intensity64_lst_lowres = zeros(1, 14);
intensity97_lst_lowres = zeros(1, 14);

intensity31_lst_lowres_err = zeros(1, 14);
intensity42_lst_lowres_err = zeros(1, 14);
intensity53_lst_lowres_err = zeros(1, 14);
intensity64_lst_lowres_err = zeros(1, 14);
intensity97_lst_lowres_err = zeros(1, 14);

h = 0;
for i = 1:14
    temp_int31 = zeros(1, length(H_LST));
    temp_int42 = zeros(1, length(H_LST));
    temp_int53 = zeros(1, length(H_LST));
    temp_int64 = zeros(1, length(H_LST));
    temp_int97 = zeros(1, length(K_LST));
    k = 0;
    for j = 1:length(H_LST)
        if str2double(H_LST{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int31(k) = intensity31(j);
            temp_int42(k) = intensity42(j);
            temp_int53(k) = intensity53(j);
            temp_int64(k) = intensity64(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    
    k = 0;
    for j = 1:length(K_LST)
        if str2double(K_LST{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int97(k) = intensity97(j);
        end
    end
    temp_int97 = temp_int97(1:k);

    intensity31_lst_lowres(i) = mean(temp_int31);
    intensity42_lst_lowres(i) = mean(temp_int42);
    intensity53_lst_lowres(i) = mean(temp_int53);
    intensity64_lst_lowres(i) = mean(temp_int64);
    intensity97_lst_lowres(i) = mean(temp_int97);
    intensity31_lst_lowres_err(i) = std(temp_int31)/sqrt(length(temp_int31));
    intensity42_lst_lowres_err(i) = std(temp_int42)/sqrt(length(temp_int42));
    intensity53_lst_lowres_err(i) = std(temp_int53)/sqrt(length(temp_int53));
    intensity64_lst_lowres_err(i) = std(temp_int64)/sqrt(length(temp_int64));
    intensity97_lst_lowres_err(i) = std(temp_int97)/sqrt(length(temp_int97));
end

%% highres
% We have LSTs from 18 to 23 and 0 to 5
H_lst_hours = zeros(1, length(H_LST));
K_lst_hours = zeros(1, length(K_LST));
lst_hours_highres = [17:0.5:23.5 0:0.5:6.5];
intensity31_lst_highres = zeros(1, 28);
intensity42_lst_highres = zeros(1, 28);
intensity53_lst_highres = zeros(1, 28);
intensity64_lst_highres = zeros(1, 28);
intensity97_lst_highres = zeros(1, 28);

intensity31_lst_highres_err = zeros(1, 28);
intensity42_lst_highres_err = zeros(1, 28);
intensity53_lst_highres_err = zeros(1, 28);
intensity64_lst_highres_err = zeros(1, 28);
intensity97_lst_highres_err = zeros(1, 28);

h = 0;
for i = 1:28
    temp_int31 = zeros(1, length(H_LST));
    temp_int42 = zeros(1, length(H_LST));
    temp_int53 = zeros(1, length(H_LST));
    temp_int64 = zeros(1, length(H_LST));
    temp_int97 = zeros(1, length(K_LST));
    k = 0;
    for j = 1:length(H_LST)
        if str2double(H_LST{j}(12:13)) + str2double(H_LST{j}(15:16))/60 >= lst_hours_highres(i)...
                && str2double(H_LST{j}(12:13)) + str2double(H_LST{j}(15:16))/60 <= lst_hours_highres(i) + 0.5
            k = k + 1;
            temp_int31(k) = intensity31(j);
            temp_int42(k) = intensity42(j);
            temp_int53(k) = intensity53(j);
            temp_int64(k) = intensity64(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    
    k = 0;
    for j = 1:length(K_LST)
        if str2double(K_LST{j}(12:13)) + str2double(K_LST{j}(15:16))/60 >= lst_hours_highres(i)...
                && str2double(K_LST{j}(12:13)) + str2double(K_LST{j}(15:16))/60 <= lst_hours_highres(i) + 0.5
            k = k + 1;
            temp_int97(k) = intensity97(j);
        end
    end
    temp_int97 = temp_int97(1:k);

    intensity31_lst_highres(i) = mean(temp_int31);
    intensity42_lst_highres(i) = mean(temp_int42);
    intensity53_lst_highres(i) = mean(temp_int53);
    intensity64_lst_highres(i) = mean(temp_int64);
    intensity97_lst_highres(i) = mean(temp_int97);
    intensity31_lst_highres_err(i) = std(temp_int31)/sqrt(length(temp_int31));
    intensity42_lst_highres_err(i) = std(temp_int42)/sqrt(length(temp_int42));
    intensity53_lst_highres_err(i) = std(temp_int53)/sqrt(length(temp_int53));
    intensity64_lst_highres_err(i) = std(temp_int64)/sqrt(length(temp_int64));
    intensity97_lst_highres_err(i) = std(temp_int97)/sqrt(length(temp_int97));
end
%% LST summer/winter lowres
lst_hours_lowres = [17 18 19 20 21 22 23 0 1 2 3 4 5 6];
intensity31_lst_winter_lowres = zeros(1, 14);
intensity42_lst_winter_lowres = zeros(1, 14);
intensity53_lst_winter_lowres = zeros(1, 14);
intensity64_lst_winter_lowres = zeros(1, 14);
intensity97_lst_winter_lowres = zeros(1, 14);

intensity31_lst_winter_lowres_err = zeros(1, 14);
intensity42_lst_winter_lowres_err = zeros(1, 14);
intensity53_lst_winter_lowres_err = zeros(1, 14);
intensity64_lst_winter_lowres_err = zeros(1, 14);
intensity97_lst_winter_lowres_err = zeros(1, 14);

for i = 1:14
    temp_int31 = zeros(1, length(H_LST_winter));
    temp_int42 = zeros(1, length(H_LST_winter));
    temp_int53 = zeros(1, length(H_LST_winter));
    temp_int64 = zeros(1, length(H_LST_winter));
    temp_int97 = zeros(1, length(K_LST_winter));
    k = 0;
    for j = 1:length(H_LST_winter)
        if str2double(H_LST_winter{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int31(k) = intensity31_winter(j);
            temp_int42(k) = intensity42_winter(j);
            temp_int53(k) = intensity53_winter(j);
            temp_int64(k) = intensity64_winter(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    
    k = 0;
    for j = 1:length(K_LST_winter)
        if str2double(K_LST_winter{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int97(k) = intensity97_winter(j);
        end
    end
    temp_int97 = temp_int97(1:k);

    intensity31_lst_winter_lowres(i) = mean(temp_int31);
    intensity42_lst_winter_lowres(i) = mean(temp_int42);
    intensity53_lst_winter_lowres(i) = mean(temp_int53);
    intensity64_lst_winter_lowres(i) = mean(temp_int64);
    intensity97_lst_winter_lowres(i) = mean(temp_int97);
    intensity31_lst_winter_lowres_err(i) = std(temp_int31)/sqrt(length(temp_int31));
    intensity42_lst_winter_lowres_err(i) = std(temp_int42)/sqrt(length(temp_int42));
    intensity53_lst_winter_lowres_err(i) = std(temp_int53)/sqrt(length(temp_int53));
    intensity64_lst_winter_lowres_err(i) = std(temp_int64)/sqrt(length(temp_int64));
    intensity97_lst_winter_lowres_err(i) = std(temp_int97)/sqrt(length(temp_int97));
end

lst_hours_lowres = [17 18 19 20 21 22 23 0 1 2 3 4 5 6];
intensity31_lst_summer_lowres = zeros(1, 14);
intensity42_lst_summer_lowres = zeros(1, 14);
intensity53_lst_summer_lowres = zeros(1, 14);
intensity64_lst_summer_lowres = zeros(1, 14);
intensity97_lst_summer_lowres = zeros(1, 14);

intensity31_lst_summer_lowres_err = zeros(1, 14);
intensity42_lst_summer_lowres_err = zeros(1, 14);
intensity53_lst_summer_lowres_err = zeros(1, 14);
intensity64_lst_summer_lowres_err = zeros(1, 14);
intensity97_lst_summer_lowres_err = zeros(1, 14);

for i = 1:14
    temp_int31 = zeros(1, length(H_LST_summer));
    temp_int42 = zeros(1, length(H_LST_summer));
    temp_int53 = zeros(1, length(H_LST_summer));
    temp_int64 = zeros(1, length(H_LST_summer));
    temp_int97 = zeros(1, length(K_LST_summer));
    k = 0;
    for j = 1:length(H_LST_summer)
        if str2double(H_LST_summer{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int31(k) = intensity31_summer(j);
            temp_int42(k) = intensity42_summer(j);
            temp_int53(k) = intensity53_summer(j);
            temp_int64(k) = intensity64_summer(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    
    k = 0;
    for j = 1:length(K_LST_summer)
        if str2double(K_LST_summer{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int97(k) = intensity97_summer(j);
        end
    end
    temp_int97 = temp_int97(1:k);

    intensity31_lst_summer_lowres(i) = mean(temp_int31);
    intensity42_lst_summer_lowres(i) = mean(temp_int42);
    intensity53_lst_summer_lowres(i) = mean(temp_int53);
    intensity64_lst_summer_lowres(i) = mean(temp_int64);
    intensity97_lst_summer_lowres(i) = mean(temp_int97);
    intensity31_lst_summer_lowres_err(i) = std(temp_int31)/sqrt(length(temp_int31));
    intensity42_lst_summer_lowres_err(i) = std(temp_int42)/sqrt(length(temp_int42));
    intensity53_lst_summer_lowres_err(i) = std(temp_int53)/sqrt(length(temp_int53));
    intensity64_lst_summer_lowres_err(i) = std(temp_int64)/sqrt(length(temp_int64));
    intensity97_lst_summer_lowres_err(i) = std(temp_int97)/sqrt(length(temp_int97));
end


%% Local solar time eyeball
%% Monthly
eyeball97_lst_lowres = zeros(14,13);
eyeball97_lst_highres = zeros(28,13);
lst_hours_lowres = [17 18 19 20 21 22 23 0 1 2 3 4 5 6];
lst_hours_highres = [17:0.5:23.5 0:0.5:6.5];

for m = 1:13
    intensity97_m = intensity97(month([K_time{:}]) == m);
    for i = 1:14
        temp_int97 = zeros(1, length(K_LST));
        k = 0;
        for j = 1:length(intensity97)
            if month(K_time{j}) == m && str2double(K_LST{j}(12:13)) == lst_hours_lowres(i)
                k = k + 1;
                temp_int97(k) = intensity97(j);
            end
        end
        temp_int97 = temp_int97(1:k);
        eyeball97_lst_lowres(i, m) = mean(temp_int97);
    end
    
    for i = 1:28
        temp_int97 = zeros(1, length(intensity97));
        k = 0;
        for j = 1:length(intensity97)
            if month(K_time{j}) == m && str2double(K_LST{j}(12:13)) + str2double(K_LST{j}(15:16))/60 >= lst_hours_highres(i)...
                    && str2double(K_LST{j}(12:13)) + str2double(K_LST{j}(15:16))/60 <= lst_hours_highres(i) + 0.5
                k = k + 1;
                temp_int97(k) = intensity97(j);
            end
        end
        temp_int97 = temp_int97(1:k);
        eyeball97_lst_highres(i, m) = mean(temp_int97);
    end
end

%% Interpolate
interp_method = 2; % 2 or 4 (or 0,2
eyeball97_lst_lowres_interp_full = inpaint_nans(repmat(eyeball97_lst_lowres([1:14],[1:12]),3,3), interp_method);
eyeball97_lst_lowres_interp = eyeball97_lst_lowres_interp_full([15:29], [13:25]);

eyeball97_lst_highres_interp_full = inpaint_nans(repmat(eyeball97_lst_highres([1:28],[1:12]),3,3), interp_method);
eyeball97_lst_highres_interp = eyeball97_lst_highres_interp_full([28:56], [13:25]);

%% Cut off airmass > 2 and see if anything changes
% K filter
intensity97_amlt2 = intensity97(K_airmass <= 2);
relative_error97_amlt2 = relative_error97(K_airmass <= 2);
K_time_amlt2 = K_time(K_airmass <= 2);

intensity97_seasonMean_amlt2 = zeros(1,12);
intensity97_seasonwMean_amlt2 = zeros(1,12);
intensity97_seasonMean_amlt2_err = zeros(1,12);

intensity97_amgt2 = intensity97(K_airmass >= 2);
relative_error97_amgt2 = relative_error97(K_airmass >= 2);
K_time_amgt2 = K_time(K_airmass >= 2);

intensity97_seasonMean_amgt2 = zeros(1,12);
intensity97_seasonwMean_amgt2 = zeros(1,12);
intensity97_seasonMean_amgt2_err = zeros(1,12);

% H filter
intensity31_amlt2 = intensity31(H_airmass <= 2);
intensity42_amlt2 = intensity42(H_airmass <= 2);
intensity53_amlt2 = intensity53(H_airmass <= 2);
intensity64_amlt2 = intensity64(H_airmass <= 2);
relative_error31_amlt2 = relative_error31(H_airmass <= 2);
relative_error42_amlt2 = relative_error42(H_airmass <= 2);
relative_error53_amlt2 = relative_error53(H_airmass <= 2);
relative_error64_amlt2 = relative_error64(H_airmass <= 2);
H_time_amlt2 = H_time(H_airmass <= 2);

intensity31_seasonMean_amlt2 = zeros(1,12);
intensity42_seasonMean_amlt2 = zeros(1,12);
intensity53_seasonMean_amlt2 = zeros(1,12);
intensity64_seasonMean_amlt2 = zeros(1,12);
intensity31_seasonwMean_amlt2 = zeros(1,12);
intensity42_seasonwMean_amlt2 = zeros(1,12);
intensity53_seasonwMean_amlt2 = zeros(1,12);
intensity64_seasonwMean_amlt2 = zeros(1,12);
intensity31_seasonMean_amlt2_err = zeros(1,12);
intensity42_seasonMean_amlt2_err = zeros(1,12);
intensity53_seasonMean_amlt2_err = zeros(1,12);
intensity64_seasonMean_amlt2_err = zeros(1,12);

intensity31_amgt2 = intensity31(H_airmass >= 2);
intensity42_amgt2 = intensity42(H_airmass >= 2);
intensity53_amgt2 = intensity53(H_airmass >= 2);
intensity64_amgt2 = intensity64(H_airmass >= 2);
relative_error31_amgt2 = relative_error31(H_airmass >= 2);
relative_error42_amgt2 = relative_error42(H_airmass >= 2);
relative_error53_amgt2 = relative_error53(H_airmass >= 2);
relative_error64_amgt2 = relative_error64(H_airmass >= 2);
H_time_amgt2 = H_time(H_airmass >= 2);

intensity31_seasonMean_amgt2 = zeros(1,12);
intensity42_seasonMean_amgt2 = zeros(1,12);
intensity53_seasonMean_amgt2 = zeros(1,12);
intensity64_seasonMean_amgt2 = zeros(1,12);
intensity31_seasonwMean_amgt2 = zeros(1,12);
intensity42_seasonwMean_amgt2 = zeros(1,12);
intensity53_seasonwMean_amgt2 = zeros(1,12);
intensity64_seasonwMean_amgt2 = zeros(1,12);
intensity31_seasonMean_amgt2_err = zeros(1,12);
intensity42_seasonMean_amgt2_err = zeros(1,12);
intensity53_seasonMean_amgt2_err = zeros(1,12);
intensity64_seasonMean_amgt2_err = zeros(1,12);

% Lower than 2
for i = 1:12
    k = 0;
    temp_int97 = zeros(1, length(intensity97_amlt2));
    temp_relerr97 = zeros(1, length(intensity97_amlt2));
    for j = 1:length(intensity97_amlt2)
        if month(K_time_amlt2{j}) == i
            k = k + 1;
            temp_int97(k) = intensity97_amlt2(j);
            temp_relerr97(k) = relative_error97_amlt2(j);
        end
    end
    temp_int97 = temp_int97(1:k);
    temp_relerr97 = temp_relerr97(1:k);
    intensity97_seasonMean_amlt2(i) = mean(temp_int97);
    intensity97_seasonwMean_amlt2(i) = sum([temp_int97.*(1./(temp_relerr97.^2))])...
        /sum(1./temp_relerr97.^2);
    intensity97_seasonMean_amlt2_err(i) = std(temp_int97)/length(temp_int97);

    % H filter
    k = 0;
    temp_int31 = zeros(1, length(intensity31_amlt2));
    temp_int42 = zeros(1, length(intensity31_amlt2));
    temp_int53 = zeros(1, length(intensity31_amlt2));
    temp_int64 = zeros(1, length(intensity31_amlt2));
    temp_relerr31 = zeros(1, length(intensity31_amlt2));
    temp_relerr42 = zeros(1, length(intensity31_amlt2));
    temp_relerr53 = zeros(1, length(intensity31_amlt2));
    temp_relerr64 = zeros(1, length(intensity31_amlt2));
    for j = 1:length(intensity31_amlt2)
        if month(H_time_amlt2{j}) == i
            k = k + 1;
            temp_int31(k) = intensity31_amlt2(j);
            temp_int42(k) = intensity42_amlt2(j);
            temp_int53(k) = intensity53_amlt2(j);
            temp_int64(k) = intensity64_amlt2(j);
            temp_relerr31(k) = relative_error31_amlt2(j);
            temp_relerr42(k) = relative_error42_amlt2(j);
            temp_relerr53(k) = relative_error53_amlt2(j);
            temp_relerr64(k) = relative_error64_amlt2(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    temp_relerr31 = temp_relerr31(1:k);
    temp_relerr42 = temp_relerr42(1:k);
    temp_relerr53 = temp_relerr53(1:k);
    temp_relerr64 = temp_relerr64(1:k);
    intensity31_seasonMean_amlt2(i) = mean(temp_int31);
    intensity42_seasonMean_amlt2(i) = mean(temp_int42);
    intensity53_seasonMean_amlt2(i) = mean(temp_int53);
    intensity64_seasonMean_amlt2(i) = mean(temp_int64);
    intensity31_seasonwMean_amlt2(i) = sum([temp_int31.*(1./(temp_relerr31.^2))])...
        /sum(1./temp_relerr31.^2);
    intensity42_seasonwMean_amlt2(i) = sum([temp_int42.*(1./(temp_relerr42.^2))])...
        /sum(1./temp_relerr42.^2);
    intensity53_seasonwMean_amlt2(i) = sum([temp_int53.*(1./(temp_relerr53.^2))])...
        /sum(1./temp_relerr53.^2);
    intensity64_seasonwMean_amlt2(i) = sum([temp_int64.*(1./(temp_relerr64.^2))])...
        /sum(1./temp_relerr64.^2);
    intensity31_seasonMean_amlt2_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_seasonMean_amlt2_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_seasonMean_amlt2_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_seasonMean_amlt2_err(i) = std(temp_int64)/length(temp_int64);
end

% Greater than 2
for i = 1:12
    k = 0;
    temp_int97 = zeros(1, length(intensity97_amgt2));
    temp_relerr97 = zeros(1, length(intensity97_amgt2));
    for j = 1:length(intensity97_amgt2)
        if month(K_time_amgt2{j}) == i
            k = k + 1;
            temp_int97(k) = intensity97_amgt2(j);
            temp_relerr97(k) = relative_error97_amgt2(j);
        end
    end
    temp_int97 = temp_int97(1:k);
    temp_relerr97 = temp_relerr97(1:k);
    intensity97_seasonMean_amgt2(i) = mean(temp_int97);
    intensity97_seasonwMean_amgt2(i) = sum([temp_int97.*(1./(temp_relerr97.^2))])...
        /sum(1./temp_relerr97.^2);
    intensity97_seasonMean_amgt2_err(i) = std(temp_int97)/length(temp_int97);

    % H filter
    k = 0;
    temp_int31 = zeros(1, length(intensity31_amgt2));
    temp_int42 = zeros(1, length(intensity31_amgt2));
    temp_int53 = zeros(1, length(intensity31_amgt2));
    temp_int64 = zeros(1, length(intensity31_amgt2));
    temp_relerr31 = zeros(1, length(intensity31_amgt2));
    temp_relerr42 = zeros(1, length(intensity31_amgt2));
    temp_relerr53 = zeros(1, length(intensity31_amgt2));
    temp_relerr64 = zeros(1, length(intensity31_amgt2));
    for j = 1:length(intensity31_amgt2)
        if month(H_time_amgt2{j}) == i
            k = k + 1;
            temp_int31(k) = intensity31_amgt2(j);
            temp_int42(k) = intensity42_amgt2(j);
            temp_int53(k) = intensity53_amgt2(j);
            temp_int64(k) = intensity64_amgt2(j);
            temp_relerr31(k) = relative_error31_amgt2(j);
            temp_relerr42(k) = relative_error42_amgt2(j);
            temp_relerr53(k) = relative_error53_amgt2(j);
            temp_relerr64(k) = relative_error64_amgt2(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    temp_relerr31 = temp_relerr31(1:k);
    temp_relerr42 = temp_relerr42(1:k);
    temp_relerr53 = temp_relerr53(1:k);
    temp_relerr64 = temp_relerr64(1:k);
    intensity31_seasonMean_amgt2(i) = mean(temp_int31);
    intensity42_seasonMean_amgt2(i) = mean(temp_int42);
    intensity53_seasonMean_amgt2(i) = mean(temp_int53);
    intensity64_seasonMean_amgt2(i) = mean(temp_int64);
    intensity31_seasonwMean_amgt2(i) = sum([temp_int31.*(1./(temp_relerr31.^2))])...
        /sum(1./temp_relerr31.^2);
    intensity42_seasonwMean_amgt2(i) = sum([temp_int42.*(1./(temp_relerr42.^2))])...
        /sum(1./temp_relerr42.^2);
    intensity53_seasonwMean_amgt2(i) = sum([temp_int53.*(1./(temp_relerr53.^2))])...
        /sum(1./temp_relerr53.^2);
    intensity64_seasonwMean_amgt2(i) = sum([temp_int64.*(1./(temp_relerr64.^2))])...
        /sum(1./temp_relerr64.^2);
    intensity31_seasonMean_amgt2_err(i) = std(temp_int31)/length(temp_int31);
    intensity42_seasonMean_amgt2_err(i) = std(temp_int42)/length(temp_int42);
    intensity53_seasonMean_amgt2_err(i) = std(temp_int53)/length(temp_int53);
    intensity64_seasonMean_amgt2_err(i) = std(temp_int64)/length(temp_int64);
end

amlt2_data = [intensity31_seasonMean_amlt2, intensity31_seasonwMean_amlt2,...
           intensity42_seasonMean_amlt2, intensity42_seasonwMean_amlt2,intensity53_seasonMean_amlt2, intensity53_seasonwMean_amlt2,...
           intensity64_seasonMean_amlt2, intensity64_seasonwMean_amlt2, intensity97_seasonMean_amlt2, intensity97_seasonwMean_amlt2];
       
amgt2_data = [intensity31_seasonMean_amgt2, intensity31_seasonwMean_amgt2,...
           intensity42_seasonMean_amgt2, intensity42_seasonwMean_amgt2,intensity53_seasonMean_amgt2, intensity53_seasonwMean_amgt2,...
           intensity64_seasonMean_amgt2, intensity64_seasonwMean_amgt2, intensity97_seasonMean_amgt2, intensity97_seasonwMean_amgt2];

%% Remove seasonal cycle
% intensity31_noy_corr = zeros(1,365);
% 
% for i = 1:365
%     if i <= 31 
%         temp_month = 1;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 32 && i <= 59
%         temp_month = 2;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 60 && i <= 90
%         temp_month = 3;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 91 && i <= 120
%         temp_month = 4;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 121 && i <= 151
%         temp_month = 5;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 152 && i <= 181
%         temp_month = 6;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 182 && i <= 211
%         temp_month = 7;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     elseif i >= 182 && i <= 211
%         temp_month = 7;
%         intensity31_noy_corr(i) = intensity31_year(temp_month);
%     end
% end
%% Figures
%% HK
if runHfilter && runKfilter
    %% Hours after sunset
    figure;
    hpss = [0:13];
    errorbar(hpss, intensity31_hpss_lowres, intensity31_hpss_lowres_err)
    hold on;
    errorbar(hpss, intensity42_hpss_lowres, intensity42_hpss_lowres_err)
    errorbar(hpss, intensity53_hpss_lowres, intensity53_hpss_lowres_err)
    errorbar(hpss, intensity64_hpss_lowres, intensity64_hpss_lowres_err)
    errorbar(hpss, intensity97_hpss_lowres, intensity97_hpss_lowres_err)
    hold off;
    legend('31','42','53','64','97')
    xlabel('Hours past sunset')
    ylabel('Intensity [a.u.]')
    title('Low resolution')
    
figure;
hold on;
box on;
    errorbar(1:12,intensity31_year_weighted, err31_year, '.-')
    errorbar(1:12,intensity42_year_weighted, err42_year, '.-') 
    errorbar(1:12,intensity53_year_weighted, err53_year, '.-') 
    errorbar(1:12,intensity64_year_weighted, err64_year, '.-')
    errorbar(1:12,intensity97_year_weighted, err97_year, '.-')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    xticks([1:12])
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt','Nov','Dec'})
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$', 'Q$(9,7)$')
    title('Weighted mean')
hold off;

figure;
hold on;
box on;
    errorbar(1:12,intensity31_year, err31_year, '.-')
    errorbar(1:12,intensity42_year, err42_year, '.-') 
    errorbar(1:12,intensity53_year, err53_year, '.-') 
    errorbar(1:12,intensity64_year, err64_year, '.-')
    errorbar(1:12,intensity97_year, err97_year, '.-')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    xticks([1:12])
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt','Nov','Dec'})
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$', 'Q$(9,7)$')
    title('Arithmetic mean')
hold off;

figure;
hold on;
box on;
    errorbar(1:12,intensity31_year./mean(intensity31_year, 'omitnan') -1, err31_year./mean(intensity31_year, 'omitnan'), '.-')
    errorbar(1:12,intensity42_year./mean(intensity42_year, 'omitnan') -1, err42_year./mean(intensity42_year, 'omitnan'), '.-') 
    errorbar(1:12,intensity53_year./mean(intensity53_year, 'omitnan') -1, err53_year./mean(intensity53_year, 'omitnan'), '.-')
    errorbar(1:12,intensity64_year./mean(intensity64_year, 'omitnan') -1, err64_year./mean(intensity64_year, 'omitnan'), '.-')
    errorbar(1:12,intensity97_year./mean(intensity97_year, 'omitnan') -1, err97_year./mean(intensity97_year, 'omitnan'), '.-')
    xlabel('Month')
    ylabel('Relative deviation from annual mean')
    xticks(1:12)
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt','Nov','Dec'})
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$', 'Q$(9,7)$')
    title('NOTCam')
    yticks([-0.6:0.1:0.6])
    ylim([-0.6, 0.6])
hold off;
    
figure;
    hpss = [0:0.5:13.5];
    errorbar(hpss, intensity31_hpss_highres, intensity31_hpss_highres_err)
    hold on;
    errorbar(hpss, intensity42_hpss_highres, intensity42_hpss_highres_err)
    errorbar(hpss, intensity53_hpss_highres, intensity53_hpss_highres_err)
    errorbar(hpss, intensity64_hpss_highres, intensity64_hpss_highres_err)
    errorbar(hpss, intensity97_hpss_highres, intensity97_hpss_highres_err)
    hold off;
    legend('31','42','53','64','97')
    xlabel('Hours past sunset')
    ylabel('Intensity [a.u.]')
    title('High resolution')
    
    %% solar depression angles
    highresdegs = 1.25:2.5:88.75;
    lowresdegs = 2.5:5:87.5;
    
    %% Low resolution | 5 deg
    figure;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity31_sdep_setting_lowres fliplr(intensity31_sdep_rising_lowres)],...
        [intensity31_sdep_setting_lowres_err fliplr(intensity31_sdep_rising_lowres_err)])
    hold on;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity42_sdep_setting_lowres fliplr(intensity42_sdep_rising_lowres)],...
        [intensity42_sdep_setting_lowres_err fliplr(intensity42_sdep_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity53_sdep_setting_lowres fliplr(intensity53_sdep_rising_lowres)],...
        [intensity53_sdep_setting_lowres_err fliplr(intensity53_sdep_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity64_sdep_setting_lowres fliplr(intensity64_sdep_rising_lowres)],...
        [intensity64_sdep_setting_lowres_err fliplr(intensity64_sdep_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity97_sdep_setting_lowres fliplr(intensity97_sdep_rising_lowres)],...
        [intensity97_sdep_setting_lowres_err fliplr(intensity97_sdep_rising_lowres_err)])
    hold off;
    xticks([0:10:180])
    xticklabels([0:10:90 80:-10:0])
    title('Combining seasons. Lowres. Midnight = 6 hpss')
    legend('31','42','53', '64', '97')
    xlabel('Solar depression angle [degrees]')
    ylabel('Intensity [a.u.]')
    
    %% Low resolution | Winter - Summer
    subplot(1,2,1)
    %figure;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity31_sdep_winter_setting_lowres fliplr(intensity31_sdep_winter_rising_lowres)],...
        [intensity31_sdep_winter_setting_lowres_err fliplr(intensity31_sdep_winter_rising_lowres_err)])
    hold on;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity42_sdep_winter_setting_lowres fliplr(intensity42_sdep_winter_rising_lowres)],...
        [intensity42_sdep_winter_setting_lowres_err fliplr(intensity42_sdep_winter_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity53_sdep_winter_setting_lowres fliplr(intensity53_sdep_winter_rising_lowres)],...
        [intensity53_sdep_winter_setting_lowres_err fliplr(intensity53_sdep_winter_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity64_sdep_winter_setting_lowres fliplr(intensity64_sdep_winter_rising_lowres)],...
        [intensity64_sdep_winter_setting_lowres_err fliplr(intensity64_sdep_winter_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity97_sdep_winter_setting_lowres fliplr(intensity97_sdep_winter_rising_lowres)],...
        [intensity97_sdep_winter_setting_lowres_err fliplr(intensity97_sdep_winter_rising_lowres_err)])
    hold off;
    xticks([0:10:180])
    xticklabels([0:10:90 80:-10:0])
    title('Winter season. Lowres. Midnight = 6.5 hpss')
    legend('31','42','53', '64', '97')
    xlabel('Solar depression angle [degrees]')
    ylabel('Intensity [a.u.]')
    
    subplot(1,2,2)
    %figure;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity31_sdep_summer_setting_lowres fliplr(intensity31_sdep_summer_rising_lowres)],...
        [intensity31_sdep_summer_setting_lowres_err fliplr(intensity31_sdep_summer_rising_lowres_err)])
    hold on;
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity42_sdep_summer_setting_lowres fliplr(intensity42_sdep_summer_rising_lowres)],...
        [intensity42_sdep_summer_setting_lowres_err fliplr(intensity42_sdep_summer_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity53_sdep_summer_setting_lowres fliplr(intensity53_sdep_summer_rising_lowres)],...
        [intensity53_sdep_summer_setting_lowres_err fliplr(intensity53_sdep_summer_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity64_sdep_summer_setting_lowres fliplr(intensity64_sdep_summer_rising_lowres)],...
        [intensity64_sdep_summer_setting_lowres_err fliplr(intensity64_sdep_summer_rising_lowres_err)])
    errorbar([2.5:5:90 (2.5:5:90)+90], [intensity97_sdep_summer_setting_lowres fliplr(intensity97_sdep_summer_rising_lowres)],...
        [intensity97_sdep_summer_setting_lowres_err fliplr(intensity97_sdep_summer_rising_lowres_err)])
    hold off;
    xticks([0:10:180])
    xticklabels([0:10:90 80:-10:0])
    title('Summer season. Lowres. Midnight = 5.5 hpss')
    legend('31','42','53', '64', '97')
    xlabel('Solar depression angle [degrees]')
    ylabel('Intensity [a.u.]')

    %% High resolution | 2.5 deg
    figure;
        plot(1.25:2.5:88.75, intensity31_sdep_setting, '-', 'Color', [0, 0.4470, 0.7410])
        hold on;
            plot(1.25:2.5:88.75, intensity31_sdep_rising, '--', 'Color', [0, 0.4470, 0.7410])
            plot(1.25:2.5:88.75, intensity42_sdep_setting, '-', 'Color', [0.8500, 0.3250, 0.0980])
            plot(1.25:2.5:88.75, intensity42_sdep_rising, '--', 'Color', [0.8500, 0.3250, 0.0980])
            plot(1.25:2.5:88.75, intensity97_sdep_setting, '-', 'Color', [0.4660, 0.6740, 0.1880])
            plot(1.25:2.5:88.75, intensity97_sdep_rising, '--', 'Color', [0.4660, 0.6740, 0.1880])
        xlabel('Solar depression angle [degrees]')
        ylabel('intensity [a.u.]')
        title('high resolution')
        legend('Q$(3,1)$ early', 'Q$(3,1)$ late','Q$(4,2)$ early', 'Q$(4,2)$ late','Q$(9,7)$ early', 'Q$(9,7)$ late')
        hold off;
    
    %% hpss vs sdep
    figure;
        plot(hPostSunSet, H_sdep, '.', K_hPostSunSet, K_sdep, '.')
        xlabel('Hours after sunset [hours]')
        ylabel('Solar depression angle [degrees]')
        legend('H band', 'K band')
        title('Solar depression angle vs. hours after sunset at NOT')
    
    %% Airmass
    figure;
        subplot(2,1,1)
        histogram(H_airmass)
        xlabel('airmass')
        title('H band')
        subplot(2,1,2)
        histogram(K_airmass)
        xlabel('airmass')
        title('K band')
        
    %% Airmass <= 2
    figure;
        plot(intensity97_year, intensity97_seasonMean_amlt2, '.')
        hold on
        plot(intensity31_year, intensity31_seasonMean_amlt2, '.')
        plot(intensity42_year, intensity42_seasonMean_amlt2, '.')
        plot(intensity53_year, intensity53_seasonMean_amlt2, '.')
        plot(intensity64_year, intensity64_seasonMean_amlt2, '.')
        hold off
        legend('97', '31', '42', '53', '64')
        title('All data vs airmass <=2 for annual arithmetic mean')
        xlabel('Intensity using full data set')
        ylabel('Intensity using only data with airmass <= 2')
        
    figure;
       plot([intensity31_year, intensity31_year_weighted, intensity42_year, intensity42_year_weighted,...
           intensity53_year, intensity53_year_weighted, intensity64_year, intensity64_year_weighted,...
           intensity97_year, intensity97_year_weighted], [intensity31_seasonMean_amlt2, intensity31_seasonwMean_amlt2,...
           intensity42_seasonMean_amlt2, intensity42_seasonwMean_amlt2,intensity53_seasonMean_amlt2, intensity53_seasonwMean_amlt2,...
           intensity64_seasonMean_amlt2, intensity64_seasonwMean_amlt2, intensity97_seasonMean_amlt2, intensity97_seasonwMean_amlt2], '.') 
       xlabel('Monthly means and weighted means. All data.')
       ylabel('Monthly means and weighted means. Airmass $<= 2$')
       title('Full data set vs data set with airmass less than 2')
       
       
    figure;
       plot(amlt2_data, amgt2_data, '.') 
       xlabel('Monthly means and weighted means. airmass $<= 2$')
       ylabel('Monthly means and weighted means. Airmass $>= 2$')
       title('Airmass gt 2 vs lt 2')
    
    %% UTC vs hpss
    figure;
        plot([hour([time{:}]) + minute([time{:}])./60], H_sdep, '.')
        
    %% Eyeball 97 sdep
    figure;
        h = pcolor(1:13, (0:5:175), eyeball_sdep_97);%, 50, 'linecolor', 'none');
        shading flat;
        colormap('Hot(1200)')
        set(h,'alphadata',~isnan(eyeball_sdep_97))
        set(gca,'color',[0.15,0.15,0.15]) 
        xticks([1:12] + 0.5)
        yticks([0:10:180])
        yticklabels([0:10:90 80:-10:0])
        xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
        xlabel('Month')
        ylabel('Solar depression angle [degrees]')
        title('Solar depression angle. 97. Time increasing upwards.')
        colorbar;
        savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_sdep_moy.fig')
        
    %% LST
    figure;
    errorbar(intensity31_lst_lowres, intensity31_lst_lowres_err)
    hold on;
    errorbar(intensity42_lst_lowres, intensity42_lst_lowres_err)
    errorbar(intensity53_lst_lowres, intensity53_lst_lowres_err)
    errorbar(intensity64_lst_lowres, intensity64_lst_lowres_err)
    errorbar(intensity97_lst_lowres, intensity97_lst_lowres_err)
    hold off;
    xticks([1:14])
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    legend('31', '42', '53', '64', '97')
    title('Local solar time, resolved in hours')
    ylim([1, 10]*1e18)
    xlim([1,14])
    
    figure;
    errorbar(intensity31_lst_highres, intensity31_lst_highres_err)
    hold on;
    errorbar(intensity42_lst_highres, intensity42_lst_highres_err)
    errorbar(intensity53_lst_highres, intensity53_lst_highres_err)
    errorbar(intensity64_lst_highres, intensity64_lst_highres_err)
    errorbar(intensity97_lst_highres, intensity97_lst_highres_err)
    hold off;
    xticks([1:2:28])
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    legend('31', '42', '53', '64', '97')
    title('Local solar time, resolved in half-hours')
    
    figure;
    errorbar(intensity31_lst_winter_lowres, intensity31_lst_winter_lowres_err)
    hold on;
    errorbar(intensity42_lst_winter_lowres, intensity42_lst_winter_lowres_err)
    errorbar(intensity53_lst_winter_lowres, intensity53_lst_winter_lowres_err)
    errorbar(intensity64_lst_winter_lowres, intensity64_lst_winter_lowres_err)
    errorbar(intensity97_lst_winter_lowres, intensity97_lst_winter_lowres_err)
    hold off;
    xticks([1:14])
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    legend('31', '42', '53', '64', '97')
    title('Local solar time, winter')
    ylim([1, 10]*1e18)
    xlim([1,14])
    
    figure;
    errorbar(intensity31_lst_summer_lowres, intensity31_lst_summer_lowres_err)
    hold on;
    errorbar(intensity42_lst_summer_lowres, intensity42_lst_summer_lowres_err)
    errorbar(intensity53_lst_summer_lowres, intensity53_lst_summer_lowres_err)
    errorbar(intensity64_lst_summer_lowres, intensity64_lst_summer_lowres_err)
    errorbar(intensity97_lst_summer_lowres, intensity97_lst_summer_lowres_err)
    hold off;
    xticks([1:14])
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    legend('31', '42', '53', '64', '97')
    title('Local solar time, summer')
    ylim([1, 10]*1e+18)
    xlim([1,14])
    
%% LST eyeball    
%% lowres
    figure;
        h = pcolor(1:13, 1:14, eyeball97_lst_lowres);%, 50, 'linecolor', 'none');
        shading flat;
        colormap('Hot(1200)')
        xticks([1:12] + 0.5)
        yticks([1:14])
        yticklabels(lst_hours_lowres)
        xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
        xlabel('Month')
        ylabel('Local solar time')
        title('LST. 97. Time increasing upwards.')
        colorbar;
        %%
        set(h,'alphadata',~isnan(eyeball_sdep_97))
        set(gca,'color',[0.15,0.15,0.15]) 
        set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
        [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
        annotation('rectangle',...
        [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
        'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
        annotation('textbox',...
        [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
        'String',{'NaN'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FitBoxToText','off');
        savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_lst_lowres.fig')
        %%
%% lowres interp        
    figure;
    contourf(eyeball97_lst_lowres_interp, 25, 'linestyle','none')
    fig_eye_lst_lowres = gca;
    colormap('Hot(1200)')
    xticks([1:12] + 0.5)
    yticks([1:14] + 0.5)
    yticklabels(lst_hours_lowres)
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('local solar time')
    shading(shade_style);
    title('lst. lowres. 97. interp.')
    colorbar;
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_lst_lowres_interp.fig')
    
 %% highres   
    figure;
        h = pcolor(eyeball97_lst_highres);%, 50, 'linecolor', 'none');
        shading flat;
        colormap('Hot(1200)')
        xticks([1:12] + 0.5)
        yticks([1:2:28])
        yticklabels(lst_hours_lowres)
        xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
        xlabel('Month')
        ylabel('Local solar time')
        title('LST. Highres. 97. Time increasing upwards.')
        colorbar;
        %%
        set(h,'alphadata',~isnan(eyeball_sdep_97))
        set(gca,'color',[0.15,0.15,0.15]) 
        set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
        [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
        annotation('rectangle',...
        [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
        'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
        annotation('textbox',...
        [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
        'String',{'NaN'},...
        'LineStyle','none',...
        'Interpreter','latex',...
        'FitBoxToText','off');
        savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_lst_highres.fig')
        %%
%% highres  interp      
    figure;
    contourf(eyeball97_lst_highres_interp, 25, 'linestyle','none')
    fig_eye_lst_highres = gca;
    colormap('Hot(1200)')
    xticks([1:12] + 0.5)
    yticks([1:2:28] + 0.5)
    yticklabels(lst_hours_lowres)
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('local solar time')
    shading(shade_style);
    title('lst. highres. 97. interp.')
    colorbar;
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_lst_highres_interp.fig')
    

end
%% H
if runHfilter
k = length(intensity31);    
%% Whole data set vs UTC
figure;
    plot([time{1:k}],[intensity31(1:k)],'.-',[time{1:k}],[intensity42(1:k)],'.-',...
        [time{1:k}],[intensity53(1:k)],'.-',[time{1:k}],[intensity64(1:k)],'.-')
    xlabel('Date')
    ylabel('Intensity [a.u.]')
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    title('Whole data set')

%% Whole data set vs index

figure;
    plot([1:k],[intensity31(1:k)],'.-',[1:k],[intensity42(1:k)],'.-',...
        [1:k],[intensity53(1:k)],'.-',[1:k],[intensity64(1:k)],'.-')
    xlabel('Index')
    ylabel('Intensity [a.u.]')
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    title('Whole data set')

%% Whole data set with error
figure;
    errorbar(1:k,[intensity31(1:k)],abserr31(1:k), '.-')
    hold on;
    errorbar(1:k,[intensity42(1:k)],abserr42(1:k), '.-')
    errorbar(1:k,[intensity53(1:k)],abserr53(1:k), '.-')
    errorbar(1:k,[intensity64(1:k)],abserr64(1:k), '.-')
    hold off;
    %datetick('x', 'yyyy/mm/dd HH:MM:SS')
    xlabel('Index')
    ylabel('Intensity')
    legend('Q31', 'Q42','Q53','Q64')
    title('Whole data set')
%% Data vs Time with error
figure;
    errorbar(datenum([time{1:k}]),[intensity31(1:k)],abserr31(1:k), '.-')
    hold on;
    errorbar(datenum([time{1:k}]),[intensity42(1:k)],abserr42(1:k), '.-')
    errorbar(datenum([time{1:k}]),[intensity53(1:k)],abserr53(1:k), '.-')
    errorbar(datenum([time{1:k}]),[intensity64(1:k)],abserr64(1:k), '.-')
    hold off;
    datetick('x', 'yyyy/mm/dd HH:MM:SS')
    xlabel('Date and time [UTC]')
    ylabel('Intensity [a.u.]')
    legend('Q$(3,1)$', 'Q$(4,2)$','Q$(5,3)$','Q$(6,4)$')
    
%% Data vs dayofyear with error   
figure;
    errorbar(day([time{1:k}], 'dayofyear'),[intensity31(1:k)],abserr31(1:k), '.')
    hold on;
    errorbar(day([time{1:k}], 'dayofyear'),[intensity42(1:k)],abserr42(1:k), '.')
    errorbar(day([time{1:k}], 'dayofyear'),[intensity53(1:k)],abserr53(1:k), '.')
    errorbar(day([time{1:k}], 'dayofyear'),[intensity64(1:k)],abserr64(1:k), '.')
    hold off;
    %datetick('x', 'yyyy/mm/dd HH:MM:SS')
    xlabel('Day of year')
    ylabel('Intensity [a.u.]')
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
%% Distribution of doys   
figure;
    histogram2([day([H_time{:}], 'dayofyear')], [year([H_time{:}])],[0:365], [2008:2017], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year')
    ylabel('Year')
    colorbar
    yticks(2008.5:2016.5)
    ytickslabels(2008:2016)
    
figure;
    histogram2([day([K_time{:}], 'dayofyear')], [year([K_time{:}])],[0:365], [2008:2017], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year')
    ylabel('Year')
    colorbar
    yticks(2008.5:2016.5)
    yticklabels(2008:2016)
%% I/I_31 vs index    
figure;
    plot([1:k],ones(1,k),[1:k], i42,'.-',[1:k], i53,'.-',[1:k], i64,'.-')
    xlabel('Index')
    ylabel('$I/I_{(3,1)}$')
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    
%% Annual mean H filter    
figure;
    hold on;
    errorbar(1:12,intensity31_year, err31_year, '.-')
    errorbar(1:12,intensity42_year, err42_year, '.-')
    errorbar(1:12,intensity53_year, err53_year, '.-')
    errorbar(1:12,intensity64_year, err64_year, '.-')
    hold off;
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    title('Global mean')
    
%% Annual weighted mean H filter    
figure;
hold on;
box on;
    errorbar(1:12,intensity31_year_weighted, err31_year_weighted, '.-')
    errorbar(1:12,intensity42_year_weighted, err42_year_weighted, '.-') 
    errorbar(1:12,intensity53_year_weighted, err53_year_weighted, '.-') 
    errorbar(1:12,intensity64_year_weighted, err64_year_weighted, '.-')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    xticks([1:12])
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt','Nov','Dec'})
    legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    title('Weighted mean')
hold off;

%% I/I_31 vs index with error
figure;
box on;
    hold on;
    plot([1:k], [ones(1,k)], '.-')
    errorbar(1:k, i42, relerr_i42.*i42, '.-')
    errorbar(1:k, i53, relerr_i53.*i53, '.-')
    errorbar(1:k, i64, relerr_i64.*i64, '.-')
    legend('Q$(3,1$)', 'Q$(4,2)$','Q$(5,3)$','Q$(6,4)$')
    ylabel('$I/I_{(3,1)}$', 'Interpreter', 'latex')
    xlabel('Index')
    title('yolo')
    hold off;

%% Monthly distributions
figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 1)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Jan')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 2)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Feb')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 3)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Mar')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 4)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Apr')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 5)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('May')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 6)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Jun')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 7)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Jul')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 8)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Aug')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 9)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Sep')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 10)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Oct')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 11)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Nov')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(hPostSunSet(month([time{:}]) == 12)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('Dec')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)
end
%% K
if runKfilter
    %% annual mean 97
    figure;
    hold on;
    errorbar(1:12,intensity97_year, err97_year, '.-')
    hold off;
    %legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    title('Global mean 97')
    
    %% Annual weighted mean
    figure;
    hold on;
    errorbar(1:12,intensity97_year_weighted, err97_year_weighted, '.-')
    hold off;
    %legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    xlabel('Month')
    ylabel('Intensity [a.u.]')
    title('Weighted Global mean 97')
    
    %% 97 vs time
    figure;
    plot([K_time{:}],[intensity97],'.-')
    xlabel('Date')
    ylabel('Intensity [a.u.]')
    %legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    title('Whole data set 97')
  %% time with error  
    figure;
    errorbar(datenum([K_time{:}]),[intensity97],abserr97, '.-')
    datetick('x', 'yyyy/mm/dd HH:MM:SS')
    xlabel('Date and time [UTC]')
    ylabel('Intensity [a.u.]')
    legend('Q$(9,7)$')
%% 97 vs index
figure;
    plot([1:K_k],[intensity97],'.-')
    xlabel('Index')
    ylabel('Intensity [a.u.]')
    %legend('Q$(3,1)$', 'Q$(4,2)$', 'Q$(5,3)$', 'Q$(6,4)$')
    title('Whole data set 97')
%% 97 vs index with error    
figure;
    errorbar(1:K_k,[intensity97],abserr97, '.-')
    hold on;
    hold off;
    %datetick('x', 'yyyy/mm/dd HH:MM:SS')
    xlabel('Index')
    ylabel('Intensity')
    %legend('Q31', 'Q42','Q53','Q64')
    title('Whole data set 97')
    
%% Hours after sunset distribution

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 1)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Jan')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 2)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Feb')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 3)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Mar')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 4)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Apr')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 5)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: May')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 6)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Jun')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 7)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Jul')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 8)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Aug')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 9)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Sep')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 10)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Oct')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 11)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Nov')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

figure('DefaultAxesFontSize',15);
histogram(floor(K_hPostSunSet(month([K_time{:}]) == 12)), [0,1,2,3,4,5,6,7,8,9,10,11,12])
xlim([0,12])
ylim([0,30])
title('K-band: Dec')
xlabel('Hours after sunset', 'FontSize', 20)
ylabel('Number of observations', 'FontSize', 20)

%% Eye ball figures
color_style = 'hot(600)'; %'default';
shade_style = 'interp'; %'interp'; %'flat';
%% Eyeball 97 hpss lowres
figure;
    pcolor(1:13, 1:13, eyeball97_moy_lowres([2:14], [1:13]))
    fig_eye = gca;
    colormap(color_style)
    xticks([1:12] + 0.5)
    yticks([1:14] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading('flat');
    title('mean of hours after sunset, for each month')
    colorbar;
    %%
    set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
    [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
    annotation('rectangle',...
    [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
    'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
    annotation('textbox',...
    [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
    'String',{'NaN'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_moy_lowres.fig')

%% 97 eyeball lowres interp
figure;
    contourf(1:13, 1:13, eyeball97_moy_lowres_interp([2:14], [1:13]), 25, 'linestyle', 'none')
    fig_eye = gca;
    colormap(color_style)
    xticks([1:12] + 0.5)
    yticks([1:14] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading(shade_style);
    title('mean of hours after sunset, for each month. 97. interp.')
    colorbar;
    %%
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_moy_lowres_interp.fig')
%% eyeball 97 highres     
figure;
    pcolor(1:13, 3:25, eyeball97_moy_highres([3:25], [1:13]))%, 25, 'linestyle', 'none')
    fig_eye_half = gca;
    colormap(color_style)
    xticks([1:12] + 0.5)
    yticks([3:2:27] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading('flat');
    title('mean of half-hours after sunset, for each month. 97.')
    colorbar;
    %%
    set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
    [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
    annotation('rectangle',...
    [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
    'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
    annotation('textbox',...
    [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
    'String',{'NaN'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_moy_highres.fig')
%% Eyeball 97 highres interp
figure;
    contourf(1:13, 3:25, eyeball97_moy_highres_interp([3:25], [1:13]), 25, 'linestyle', 'none')
    fig_eye_half = gca;
    colormap(color_style)
    xticks([1:12] + 0.5)
    yticks([3:2:27] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading(shade_style);
    title('mean of half-hours after sunset, for each month. 97. interp')
    colorbar;
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_moy_highres_interp.fig')
%% eyeball lowres 97 noy    
figure;
    pcolor(1:366, 1:13, eyeball97_noy_lowres([2:14], [1:366]))
    fig_eye_noy = gca;
    colormap(color_style)
    xticks([1:32:365] + 0.5)
    yticks([1:14] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading('flat');
    title('mean of hours after sunset, for each day. 97.')
    colorbar;
%%
    set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
    [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
    annotation('rectangle',...
    [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
    'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
    annotation('textbox',...
    [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
    'String',{'NaN'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_noy_lowres.fig')
%% eyeball 97 lowres noy interp   
figure;
    contourf(1:366, 1:13, eyeball97_noy_lowres_interp([2:14], [1:366]), 25, 'linestyle', 'none')
    fig_eye_noy = gca;
    colormap(color_style)
    xticks([1:32:365] + 0.5)
    yticks([1:14] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading(shade_style);
    title('mean of hours after sunset, for each day.97.interp')
    colorbar;
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_noy_lowres_interp.fig')
%% eyeball 97 highres noy    
figure;
    pcolor(1:366, 3:25, eyeball97_noy_highres([3:25], [1:366]))
    fig_eye_noy_half = gca;
    colormap(color_style)
    xticks([1:32:365] + 0.5)
    yticks([3:2:27] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading('flat');
    title('mean of half-hours after sunset, for each day. 97.')
    colorbar;
%%
    set(gca,'color',[0.149019607843137 0.149019607843137 0.149019607843137], 'colormap',...
    [0.331441291040623 0 0;0.335492321823287 0 0;0.339543352605952 0 0;0.343594383388616 0 0;0.34764541417128 0 0;0.351696444953945 0 0;0.355747475736609 0 0;0.359798506519273 0 0;0.363849537301938 0 0;0.367900568084602 0 0;0.371951598867266 0 0;0.376002629649931 0 0;0.380053660432595 0 0;0.384104691215259 0 0;0.388155721997923 0 0;0.392206752780588 0 0;0.396257783563252 0 0;0.400308814345916 0 0;0.404359845128581 0 0;0.408410875911245 0 0;0.412461906693909 0 0;0.416512937476574 0 0;0.420563968259238 0 0;0.424614999041902 0 0;0.428666029824567 0 0;0.432717060607231 0 0;0.436768091389895 0 0;0.440819122172559 0 0;0.444870152955224 0 0;0.448921183737888 0 0;0.452972214520552 0 0;0.457023245303217 0 0;0.461074276085881 0 0;0.465125306868545 0 0;0.46917633765121 0 0;0.473227368433874 0 0;0.477278399216538 0 0;0.481329429999203 0 0;0.485380460781867 0 0;0.489431491564531 0 0;0.493482522347195 0 0;0.49753355312986 0 0;0.501584583912524 0 0;0.505635614695189 0 0;0.509686645477853 0 0;0.513737676260517 0 0;0.517788707043181 0 0;0.521839737825846 0 0;0.52589076860851 0 0;0.529941799391174 0 0;0.533992830173839 0 0;0.538043860956503 0 0;0.542094891739167 0 0;0.546145922521832 0 0;0.550196953304496 0 0;0.55424798408716 0 0;0.558299014869825 0 0;0.562350045652489 0 0;0.566401076435153 0 0;0.570452107217818 0 0;0.574503138000482 0 0;0.578554168783146 0 0;0.582605199565811 0 0;0.586656230348475 0 0;0.59070726113114 0 0;0.594758291913804 0 0;0.598809322696469 0 0;0.602860353479133 0 0;0.606911384261798 0 0;0.610962415044462 0 0;0.615013445827126 0 0;0.619064476609791 0 0;0.623115507392455 0 0;0.62716653817512 0 0;0.631217568957784 0 0;0.635268599740449 0 0;0.639319630523113 0 0;0.643370661305777 0 0;0.647421692088442 0 0;0.651472722871106 0 0;0.655523753653771 0 0;0.659574784436435 0 0;0.6636258152191 0 0;0.667676846001764 0 0;0.671727876784428 0 0;0.675778907567093 0 0;0.679829938349757 0 0;0.683880969132422 0 0;0.687931999915086 0 0;0.691983030697751 0 0;0.696034061480415 0 0;0.70008509226308 0 0;0.704136123045744 0 0;0.708187153828408 0 0;0.712238184611073 0 0;0.716289215393737 0 0;0.720340246176402 0 0;0.724391276959066 0 0;0.728442307741731 0 0;0.732493338524395 0 0;0.736544369307059 0 0;0.740595400089724 0 0;0.744646430872388 0 0;0.748697461655053 0 0;0.752748492437717 0 0;0.756799523220381 0 0;0.760850554003046 0 0;0.76490158478571 0 0;0.768952615568375 0 0;0.773003646351039 0 0;0.777054677133704 0 0;0.781105707916368 0 0;0.785156738699033 0 0;0.789207769481697 0 0;0.793258800264361 0 0;0.797309831047026 0 0;0.80136086182969 0 0;0.805411892612355 0 0;0.809462923395019 0 0;0.813513954177684 0 0;0.817564984960348 0 0;0.821616015743012 0 0;0.825667046525677 0 0;0.829718077308341 0 0;0.833769108091006 0 0;0.83782013887367 0 0;0.841871169656334 0 0;0.845922200438999 0 0;0.849973231221663 0 0;0.854024262004328 0 0;0.858075292786992 0 0;0.862126323569656 0 0;0.866177354352321 0 0;0.870228385134985 0 0;0.874279415917649 0 0;0.878330446700314 0 0;0.882381477482978 0 0;0.886432508265642 0 0;0.890483539048307 0 0;0.894534569830971 0 0;0.898585600613636 0 0;0.9026366313963 0 0;0.906687662178964 0 0;0.910738692961629 0 0;0.914789723744293 0 0;0.918840754526957 0 0;0.922891785309622 0 0;0.926942816092286 0 0;0.930993846874951 0 0;0.935044877657615 0 0;0.939095908440279 0 0;0.943146939222944 0 0;0.947197970005608 0 0;0.951249000788272 0 0;0.955300031570937 0 0;0.959351062353601 0 0;0.963402093136266 0 0;0.96745312391893 0 0;0.971504154701594 0 0;0.975555185484259 0 0;0.979606216266923 0 0;0.983657247049588 0 0;0.987708277832252 0 0;0.991759308614916 0 0;0.995806899951635 3.43944594599642e-06 0;0.998789593369118 0.00107177681112655 0;0.999865498621048 0.00404690234186127 0;0.99999802076514 0.00796541098043363 0;1 0.0120144625282381 0;1 0.0160654933109025 0;1 0.0201165240935669 0;1 0.0241675548762312 0;1 0.0282185856588956 0;1 0.03226961644156 0;1 0.0363206472242244 0;1 0.0403716780068887 0;1 0.0444227087895531 0;1 0.0484737395722175 0;1 0.0525247703548819 0;1 0.0565758011375462 0;1 0.0606268319202106 0;1 0.064677862702875 0;1 0.0687288934855394 0;1 0.0727799242682037 0;1 0.0768309550508681 0;1 0.0808819858335325 0;1 0.0849330166161969 0;1 0.0889840473988612 0;1 0.0930350781815256 0;1 0.09708610896419 0;1 0.101137139746854 0;1 0.105188170529519 0;1 0.109239201312183 0;1 0.113290232094847 0;1 0.117341262877512 0;1 0.121392293660176 0;1 0.125443324442841 0;1 0.129494355225505 0;1 0.133545386008169 0;1 0.137596416790834 0;1 0.141647447573498 0;1 0.145698478356162 0;1 0.149749509138827 0;1 0.153800539921491 0;1 0.157851570704155 0;1 0.16190260148682 0;1 0.165953632269484 0;1 0.170004663052149 0;1 0.174055693834813 0;1 0.178106724617477 0;1 0.182157755400142 0;1 0.186208786182806 0;1 0.19025981696547 0;1 0.194310847748135 0;1 0.198361878530799 0;1 0.202412909313464 0;1 0.206463940096128 0;1 0.210514970878792 0;1 0.214566001661457 0;1 0.218617032444121 0;1 0.222668063226785 0;1 0.22671909400945 0;1 0.230770124792114 0;1 0.234821155574779 0;1 0.238872186357443 0;1 0.242923217140107 0;1 0.246974247922772 0;1 0.251025278705436 0;1 0.2550763094881 0;1 0.259127340270765 0;1 0.263178371053429 0;1 0.267229401836094 0;1 0.271280432618758 0;1 0.275331463401422 0;1 0.279382494184087 0;1 0.283433524966751 0;1 0.287484555749415 0;1 0.29153558653208 0;1 0.295586617314744 0;1 0.299637648097408 0;1 0.303688678880073 0;1 0.307739709662737 0;1 0.311790740445401 0;1 0.315841771228065 0;1 0.31989280201073 0;1 0.323943832793394 0;1 0.327994863576058 0;1 0.332045894358722 0;1 0.336096925141387 0;1 0.340147955924051 0;1 0.344198986706715 0;1 0.348250017489379 0;1 0.352301048272044 0;1 0.356352079054708 0;1 0.360403109837372 0;1 0.364454140620036 0;1 0.368505171402701 0;1 0.372556202185365 0;1 0.376607232968029 0;1 0.380658263750694 0;1 0.384709294533358 0;1 0.388760325316022 0;1 0.392811356098686 0;1 0.396862386881351 0;1 0.400913417664015 0;1 0.404964448446679 0;1 0.409015479229343 0;1 0.413066510012008 0;1 0.417117540794672 0;1 0.421168571577336 0;1 0.42521960236 0;1 0.429270633142665 0;1 0.433321663925329 0;1 0.437372694707993 0;1 0.441423725490657 0;1 0.445474756273322 0;1 0.449525787055986 0;1 0.45357681783865 0;1 0.457627848621314 0;1 0.461678879403979 0;1 0.465729910186643 0;1 0.469780940969307 0;1 0.473831971751971 0;1 0.477883002534636 0;1 0.4819340333173 0;1 0.485985064099964 0;1 0.490036094882628 0;1 0.494087125665293 0;1 0.498138156447957 0;1 0.502189187230621 0;1 0.506240218013285 0;1 0.51029124879595 0;1 0.514342279578614 0;1 0.518393310361278 0;1 0.522444341143942 0;1 0.526495371926607 0;1 0.530546402709271 0;1 0.534597433491935 0;1 0.5386484642746 0;1 0.542699495057264 0;1 0.546750525839928 0;1 0.550801556622592 0;1 0.554852587405257 0;1 0.558903618187921 0;1 0.562954648970585 0;1 0.567005679753249 0;1 0.571056710535913 0;1 0.575107741318578 0;1 0.579158772101242 0;1 0.583209802883906 0;1 0.587260833666571 0;1 0.591311864449235 0;1 0.595362895231899 0;1 0.599413926014563 0;1 0.603464956797228 0;1 0.607515987579892 0;1 0.611567018362556 0;1 0.61561804914522 0;1 0.619669079927885 0;1 0.623720110710549 0;1 0.627771141493213 0;1 0.631822172275877 0;1 0.635873203058542 0;1 0.639924233841206 0;1 0.64397526462387 0;1 0.648026295406534 0;1 0.652077326189199 0;1 0.656128356971863 0;1 0.660179387754527 0;1 0.664230418537191 0;1 0.668281449319856 0;1 0.67233248010252 0;1 0.676383510885184 0;1 0.680434541667849 0;1 0.684485572450513 0;1 0.688536603233177 0;1 0.692587634015841 0;1 0.696638664798506 0;1 0.70068969558117 0;1 0.704740726363834 0;1 0.708791757146498 0;1 0.713023839584252 0;1 0.716653820006239 0;1 0.720283800428227 0;1 0.723913780850214 0;1 0.727543761272201 0;1 0.731173741694188 0;1 0.734803722116176 0;1 0.738433702538163 0;1 0.74206368296015 0;1 0.745693663382137 0;1 0.749323643804125 0;1 0.752953624226112 0;1 0.756583604648099 0;1 0.760213585070086 0;1 0.763843565492073 0;1 0.767473545914061 0;1 0.771103526336048 0;1 0.774733506758035 0;1 0.778363487180022 0;1 0.78199346760201 0;1 0.785623448023997 0;1 0.789253428445984 0;1 0.792883408867971 0;1 0.796513389289959 0;1 0.800143369711946 0;1 0.803773350133933 0;1 0.80740333055592 0;1 0.811033310977908 0;1 0.814663291399895 0;1 0.818293271821882 0;1 0.821923252243869 0;1 0.825553232665856 0;1 0.829183213087844 0;1 0.832813193509831 0;1 0.836443173931818 0;1 0.840073154353805 0;1 0.843703134775793 0;1 0.84733311519778 0;1 0.850963095619767 0;1 0.854593076041754 0;1 0.858223056463742 0;1 0.861853036885729 0;1 0.865483017307716 0;1 0.869112997729703 0;1 0.872742978151691 0;1 0.876372958573678 0;1 0.880002938995665 0;1 0.883632919417652 0;1 0.887262899839639 0;1 0.890892880261627 0;1 0.894522860683614 0;1 0.898152841105601 0;1 0.901782821527588 0;1 0.905412801949576 0;1 0.909042782371563 0;1 0.91267276279355 0;1 0.916302743215537 0;1 0.919932723637525 0;1 0.923562704059512 0;1 0.927192684481499 0;1 0.930822664903486 0;1 0.934452645325474 0;1 0.938082625747461 0;1 0.941712606169448 0;1 0.945342586591435 0;1 0.948972567013423 0;1 0.95260254743541 0;1 0.956232527857397 0;1 0.959862508279384 0;1 0.963492488701371 0;1 0.967122469123359 0;1 0.970752449545346 0;1 0.974382429967333 0;1 0.97801241038932 0;1 0.981642390811308 0;1 0.985272371233295 0;1 0.988902351655282 0;1 0.992532218487635 1.7038445099368e-07;1 0.996157817343127 6.74273419432007e-06;1 0.999439815419274 0.000528716252954609;1 0.999997903876187 0.00513655420056591;1 0.999999949287327 0.0105784567168368;1 1 0.0160233512808084;1 1 0.0214683219137892;1 1 0.0269132925467701;1 1 0.0323582631797509;1 1 0.0378032338127318;1 1 0.0432482044457126;1 1 0.0486931750786935;1 1 0.0541381457116744;1 1 0.0595831163446552;1 1 0.0650280869776361;1 1 0.0704730576106169;1 1 0.0759180282435978;1 1 0.0813629988765786;1 1 0.0868079695095595;1 1 0.0922529401425403;1 1 0.0976979107755212;1 1 0.103142881408502;1 1 0.108587852041483;1 1 0.114032822674464;1 1 0.119477793307445;1 1 0.124922763940425;1 1 0.130367734573406;1 1 0.135812705206387;1 1 0.141257675839368;1 1 0.146702646472349;1 1 0.15214761710533;1 1 0.157592587738311;1 1 0.163037558371291;1 1 0.168482529004272;1 1 0.173927499637253;1 1 0.179372470270234;1 1 0.184817440903215;1 1 0.190262411536196;1 1 0.195707382169177;1 1 0.201152352802157;1 1 0.206597323435138;1 1 0.212042294068119;1 1 0.2174872647011;1 1 0.222932235334081;1 1 0.228377205967062;1 1 0.233822176600043;1 1 0.239267147233023;1 1 0.244712117866004;1 1 0.250157088498985;1 1 0.255602059131966;1 1 0.261047029764947;1 1 0.266492000397928;1 1 0.271936971030909;1 1 0.277381941663889;1 1 0.28282691229687;1 1 0.288271882929851;1 1 0.293716853562832;1 1 0.299161824195813;1 1 0.304606794828794;1 1 0.310051765461775;1 1 0.315496736094755;1 1 0.320941706727736;1 1 0.326386677360717;1 1 0.331831647993698;1 1 0.337276618626679;1 1 0.34272158925966;1 1 0.348166559892641;1 1 0.353611530525621;1 1 0.359056501158602;1 1 0.364501471791583;1 1 0.369946442424564;1 1 0.375391413057545;1 1 0.380836383690526;1 1 0.386281354323507;1 1 0.391726324956487;1 1 0.397171295589468;1 1 0.402616266222449;1 1 0.40806123685543;1 1 0.413506207488411;1 1 0.418951178121391;1 1 0.424396148754372;1 1 0.429841119387353;1 1 0.435286090020334;1 1 0.440731060653315;1 1 0.446176031286296;1 1 0.451621001919276;1 1 0.457065972552257;1 1 0.462510943185238;1 1 0.467955913818219;1 1 0.4734008844512;1 1 0.478845855084181;1 1 0.484290825717162;1 1 0.489735796350143;1 1 0.495180766983124;1 1 0.500625737616106;1 1 0.506070708249087;1 1 0.511515678882068;1 1 0.516960649515049;1 1 0.52240562014803;1 1 0.527850590781012;1 1 0.533295561413993;1 1 0.538740532046974;1 1 0.544185502679955;1 1 0.549630473312936;1 1 0.555075443945917;1 1 0.560520414578899;1 1 0.56596538521188;1 1 0.571410355844861;1 1 0.576855326477842;1 1 0.582300297110823;1 1 0.587745267743804;1 1 0.593190238376786;1 1 0.598635209009767;1 1 0.604080179642748;1 1 0.609525150275729;1 1 0.61497012090871;1 1 0.620415091541691;1 1 0.625860062174673;1 1 0.631305032807654;1 1 0.636750003440635;1 1 0.642194974073616;1 1 0.647639944706597;1 1 0.653084915339579;1 1 0.658504314154768;1 1 0.664186626681104;1 1 0.67007808937091;1 1 0.675969552060716;1 1 0.681861014750521;1 1 0.687752477440327;1 1 0.693643940130133;1 1 0.699535402819939;1 1 0.705426865509745;1 1 0.711318328199551;1 1 0.717209790889357;1 1 0.723101253579163;1 1 0.728992716268969;1 1 0.734884178958775;1 1 0.740775641648581;1 1 0.746667104338387;1 1 0.752558567028193;1 1 0.758450029717999;1 1 0.764341492407805;1 1 0.770232955097611;1 1 0.776124417787416;1 1 0.782015880477222;1 1 0.787907343167028;1 1 0.793798805856834;1 1 0.79969026854664;1 1 0.805581731236446;1 1 0.811473193926252;1 1 0.817364656616058;1 1 0.823256119305864;1 1 0.82914758199567;1 1 0.835039044685476;1 1 0.840930507375282;1 1 0.846821970065088;1 1 0.852713432754894;1 1 0.8586048954447;1 1 0.864496358134506;1 1 0.870387820824312;1 1 0.876279283514118;1 1 0.882170746203923;1 1 0.888062208893729;1 1 0.893953671583535;1 1 0.899845134273341;1 1 0.905736596963147;1 1 0.911628059652953;1 1 0.917519522342759;1 1 0.923410985032565;1 1 0.929302447722371;1 1 0.935193910412177;1 1 0.941085373101983;1 1 0.946976835791789;1 1 0.952868298481595;1 1 0.958759761171401;1 1 0.964651223861207;1 1 0.970542686551013;1 1 0.976434149240818;1 1 0.982325611930624;1 1 0.98821707462043;1 1 0.994108537310233;1 1 0.999999999999974])
    annotation('rectangle',...
    [0.829166666666665 0.0458253968253974 0.0416666666666672 0.0383015873015867],...
    'FaceColor',[0.149019607843137 0.149019607843137 0.149019607843137]);
    annotation('textbox',...
    [0.869047619047619 0.022809523809524 0.0875000000000002 0.0650793650793651],...
    'String',{'NaN'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_noy_highres.fig')
%% eyeball highres 97 noy interp
figure;
    contourf(1:366, 3:25, eyeball97_noy_highres_interp([3:25], [1:366]), 25, 'linestyle','none')
    fig_eye_noy_half = gca;
    colormap(color_style)
    xticks([1:32:365] + 0.5)
    yticks([3:2:27] + 0.5)
    yticklabels({1:14})
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    xlabel('Month')
    ylabel('Hours after sunset')
    shading(shade_style);
    title('mean of half-hours after sunset, for each day. 97. interp.')
    colorbar;
    savefig('C:\Users\birke\Documents\master\resultsSummary\results_backup\97_eyeball_noy_highres_interp.fig')
end

