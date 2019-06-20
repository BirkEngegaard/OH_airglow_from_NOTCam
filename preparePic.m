function [counts31,spec_relerr31,counts42,spec_relerr42,counts53,spec_relerr53, counts64,spec_relerr64,...
    SAT, specificDark, exposureTime, obsTime_mean, airmass, rampMode]...
    = preparePic(path,file,trans,filter,dark,hot,flat,bad, figpath, S, data, showFig, useMasterDark, nightMaster)

obsTime_mean = 0;
if useMasterDark
    specificDark = 0;
elseif nightMaster
    specificDark = 2;
else
    specificDark = 1;
end

%% final picture:
%data = fitsread(strcat(path,file,'.fits'),'Image');
%S=fitsinfo(strcat(path,file,'.fits'));
k=S.PrimaryData.Keywords;
%%
%% Get observation time and date, and gain for each subdetector
[gindex, ~] = find(strcmp(S.PrimaryData.Keywords,'GAIN1'));
[oindex, ~] = find(strcmp(S.PrimaryData.Keywords,'DATE-OBS'));
obsTime = datetime(replace(S.PrimaryData.Keywords{oindex, 2}, 'T', ' '), 'Format', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC'); % Start of observation, UTC
gain = [S.PrimaryData.Keywords{[gindex:(gindex+4)],2}]; %[kw1 - kw4] <=> [Q0, Q1, Q3, Q2]
%detecorPix = {S.PrimaryData.Keywords{133,2}, S.PrimaryData.Keywords{134,2},...
%    S.PrimaryData.Keywords{135,2},S.PrimaryData.Keywords{136,2}};

%% Get exposure time and choose the observation data point to be in the middle of observation
[eindex,~] = find(strcmp(k,'EXPMODE'));
expMode = S.PrimaryData.Keywords{eindex, 2}; % Exposure mode e.g. 'frames 3.6 3'
expMode = split(expMode, ' '); % split string into array for easier handling
try
    exposureTime = str2double(expMode{2})*str2double(expMode{3}); % Total time = time of ramp * num of ramps
    obsTime_mean = obsTime + seconds(exposureTime/2);
catch
    warning('Couldn''t determine exposure time. Default value of 10.8 used')
    exposureTime = 10.8;
    obsTime_mean = obsTime;
end

rampMode = [str2double(expMode{2}), str2double(expMode{3})]; % [time of ramp, number of ramps]

% if str2double(expMode{3}) < 3 % Throw away samples with less then 3 ramps
%     counts31 = -1;
%     counts42 = -1;
%     counts53 = -1;
%     counts64 = -1;
%     SAT = -1;
%     spec_relerr31 = -1;
%     spec_relerr42 = -1;
%     spec_relerr53 = -1;
%     spec_relerr64 = -1;
%     return
% end

%% UTC time to apparent and mean solar time
[SAT,SMT] = UTC2SolarApparentTime(replace(datestr(obsTime_mean),'-', '/'), -17.8851);
%% Hours after sunset 


%% Convert data from ADUs to counts
%{'[1:512,1:512]';'[513:1024,1:512]';'[1:512,513:1024]';'[513:1024, 513:1024]'}
% data(1:512, 1:512) = data(1:512, 1:512)*gain(1);
% data(513:1024, 1:512) = data(513:1024, 1:512)*gain(2);
% data(1:512, 513:1024) = data(1:512, 513:1024)*gain(3);
% data(513:1024, 513:1024) = data(513:1024, 513:1024)*gain(4);
% 
% figure;
% contourf(data,'LineColor','none');
% pcolor(data); shading interp;
% colormap('hot(600)')
% colorbar;

%% Convert data from ADUs to ADUs/s
data = data/exposureTime;
meanGain = 2.6378; % counts/ADU

%% Data processing
data(data < -0) = -1;
%zero

%figpath = sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s_figs\\',file);

%% Save data to plot/study later
data_plot = data;
data_plot(data_plot > 500/exposureTime) = 500/exposureTime;
dark_plot = dark;
dark_plot(dark_plot > 500/exposureTime) = 500/exposureTime;
if showFig
    figure;
    contourf(data_plot,'LineColor','none');
    pcolor(data_plot); shading interp;
    colormap('hot(600)')
    colorbar;
    title(sprintf('%s: Data raw',file))
    
    figure;
    contourf(dark_plot,'LineColor','none');
    pcolor(dark_plot); shading interp;
    colormap('hot(600)')
    colorbar;
    title(sprintf('%s: dark frame',file))
    xlabel('Spatial dimension [pixel]')
    ylabel('Spectral dimension [pixel]')
end
% dlmwrite(strcat(figpath, '\\data_unproc_sat1000.txt'), data_plot) % Save data for plotting later
%%

data(abs(data) > (10000/exposureTime)) = -1; % Handle muon 10000
star = sum(data); %Handle Star
pos = find(star>(26)*1024);
if size(pos)>0
    starpos = pos(round(length(pos)/2));
    if starpos<=994 && starpos>=31
        data(:,starpos-30:starpos+30) = -1;
    elseif starpos > 994
        data(:,(starpos-30):1024) = -1;
    elseif starpose < 31
        data(:,1:(starpos+30)) = -1;
    end
end

pix = 1:1024;
if useMasterDark
    dark = dark/20;
end
counter = 0;
while sum(sum(data,2) < sum(dark,2)) > 10
    dark = dark*0.9;
    specificDark = -1; % To indicate fishyness
    counter = counter + 1;
    if counter > 100
        break;
    end
end
if showFig
    figure;
    plot(pix,sum(data,2),'r', pix, sum(dark,2), 'b')
    legend('Raw image', 'Dark frame')
    title(sprintf('%s',file))
    xlabel('Spectral dimension [pixel]')
    ylabel('Intensity [ADU/s]')
    xlim([0,1024])
end

data_hold = data;
%% Upper region
data_4 = data_hold(768:1024,:);
data_4(abs(data_4) > (100)) = -1; % Handle muon 300
data([768:1024],:) = data_4; 

%% Upper middle
data_3 = data_hold(512:767,:);
data_3(abs(data_3) > (150)) = -1; % Handle muon 300
data([512:767],:) = data_3; 

%% Lower middle
data_2 = data_hold(256:511,:);
data_2(abs(data_2) > (150)) = -1; % Handle muon 300
data([256:511],:) = data_2; 

%% Lower
data_1 = data_hold(1:255,:);
data_1(abs(data_1) > (620)) = -1; % Handle muon 300
data([1:255],:) = data_1; 

%(abs(data) > (10000/10.8)) = -1; % Handle muon 420 in 1:100 spectral area. 

% data_500_620 = data_hold(500:620,:);
% data_500_620(abs(data_500_620) > (10000/10.8)) = -1; %Handle muon 500 in 500:620 spectral area. 
% data([500:620],:) = data_500_620;

data = data - dark; % Normalise dark so that data is larger than dark
data(bad == 1) = -1;
data(hot == 0) = -1;
bsumdata = data;
bsumdata(data>=0)=0;
bsumdata(data<0)=1;
baderased = sum(bsumdata,2);

if showFig
    data_plot = data;
    data_plot(data_plot > 500/exposureTime) = 500/exposureTime;
    data_plot(data_plot<0) = 0;
    figure;
    contourf(data_plot,'LineColor','none');
    pcolor(data_plot); shading interp;
    colormap('hot(600)')
    colorbar;
    title(sprintf('%s: Filtered data',file))
    
    figure;
    data_plot(data_plot<0) = 0;
    plot(pix,sum(data_plot,2))
    legend('data - dark - cosmic - bad')
    title(sprintf('%s',file))
end


%%
% data_plot = data;
% data_plot(data_plot > 150) = 150;
% data_plot(data_plot < 0) = -1;
% dlmwrite(strcat(figpath, '\\data_preflat_sat150.txt'), data_plot) % Save data for plotting later
%%

data(data < 0) = 0;

%flat = flat./bb;
%data = data./flat;


%% 
% data_plot = data;
% data_plot(data_plot > 150) = 150;
% data_plot(data_plot < 0) = -1;
% dlmwrite(strcat(figpath, '\\data_postflat_sat150.txt'), data_plot) % Save data for plotting later
%%
%% Linearize

ldata = zeros(size(data)); %linearise
lflat = zeros(size(flat));
for i=1:1000
    for j=1:1000
        ldata(j,i)=data(floor(f(i,j)),i); %Arcane magic, do not touch! it works! just believe me!! Kind regards Christoph from 2016
        lflat(j,i)=flat(floor(f(i,j)),i);
    end  
end

ldata_plot = ldata;
ldata_plot(ldata_plot > 500/exposureTime) = 500/exposureTime;
if showFig
    figure;
    contourf(ldata_plot,'LineColor','none');
    pcolor(ldata_plot); shading interp;
    colormap('hot(600)')
    colorbar;
    title(sprintf('%s: Linearized',file))
    xlabel('Spatial dimension [pixel]')
    ylabel('Spectral dimension [pixel]')
end
[tindex, ~] = find(strcmp(k,'TELALT'));
telalt = k{tindex,2};
telalt = deg2rad(telalt);
ldata=abs(ldata.*sin(telalt));
spec = sum(ldata,2); %add up
baderased(baderased == 1024) = 10024;
spec = spec.*1024./(1024-baderased);
spec(spec < 0) = -1;
spec = flipud(spec);
sflat = mean(lflat,2);
sflat = flipud(sflat);
sigma = flipud(std(ldata,0,2));
ldata = flipud(ldata);
lflat = flipud(lflat);
airmass = 1/sin(telalt);

%% Calibration
if filter == 'K'
    if exist(strcat(path,'lambda',filter,'.txt'),'file')==2
        lambda = dlmread(strcat(path,'lambda',filter,'.txt'));
    else
        lambda = dlmread(strcat('Calibration/lambda',filter,'.txt'));
        warning('Night not yet calibrated, sample calibration used!');
    end
else
% 
% pix = transpose([1:1024]);
% Q31_1 = [(117-30):(117+30)]; 
% Q42_1 = [(385-30):(385+30)];
% Q53_1 = [(656-30):(656+30)];
% Q64_1 = [(117-30):(117+30)];

%lam = [1.50559;1.58331;1.66923;1.76515];
spec = flipud(spec); % flip back to original spectrum
pixel = linspace(1,1024,1024)';
pixRough = [896;862;822;810;785;771;590;551;540;502;367;323;307;288;277;251;239];
pix = zeros(size(pixRough));
i=1;
me = median(spec);
for p = pixRough'
    short = spec(p-10:p+10);
    y = linspace(p-10,p+10,21)';
    gauss = fittype('a*exp(-((x-m)/s)^2)+c','coeff',{'a','m','s','c'});
    f1=fit(y,short,gauss,'StartPoint',[21000, p, 1.5, me],'Lower',[100,p-8,1,me*0.1],'Upper',[40000,p+8,4,4*me]);
%     figure;
%     plot(f1,y,short);
    coeffs = coeffvalues(f1);
    pix(i)=coeffs(2);
    i=i+1;
end
lam = [1.50555;1.51870;1.52876;1.53319;1.53951;1.54316;1.59725;1.60796;1.61281;1.62347;1.66923;1.68403;1.69033;1.69548;1.70082;1.70781;1.71229];

poly3 = fittype('a*x^3+b*x^2+c*x+d','coeff',{'a','b','c','d'});
f1=fit(pix,lam,poly3,'StartPoint',[-1.6e-11,8.5e-08,3.4e-4,1.9]);
% figure;
% plot(f1,pix,lam);

coeffs = coeffvalues(f1);

p0 = coeffs(4);%1.944;
p1 = coeffs(3);%3.035e-4;
p2 = coeffs(2);%1.536e-7;
p3 = coeffs(1);%-5.705e-11;

%lambda = zeros(size(spec));
lambda = p0 + p1*pixel + p2*pixel.^2 + p3*pixel.^3;
% transpose(lambda);

residual = lam - (p0 + p1*pix + p2*pix.^2 + p3*pix.^3);
% figure;
% plot(pix,residual);
% title('Residual');
% xlabel('pixel number');
% ylabel('\Delta\Lambda (\mum)');

% figure;
% plot(lambda,spec);

lambda = flipud(lambda); % Flip so lambda increases towards to right
spec = flipud(spec);
end
%% End calibration
wavenumber = (1./lambda)*10000; % [1/cm]

%divided by flatfield
%Clean flatfield from lamp's blackbody
%bb = 1./(lambda.^5).*(1./(exp((3e14*6.626e-34)./(1.38e-23*3200.*lambda))-1));%Because lambda is known here
bb_of_lambda = @(lambda) ((2*(10^18)*3e8)./(lambda.^4)).*(1./(exp((10^6*6.626e-34*3e+8)./(1.38e-23*3200.*lambda))-1));
%bb = ((2*(10^18)*3e8)./(lambda.^4)).*(1./(exp((10^6*6.626e-34*3e+8)./(1.38e-23*3200.*lambda))-1)); % ph/m2/sr/um
%bb = bb/meanGain; % to go from counts to ADU
%bb = repmat(bb/sum(bb)*1024,1,1024);
bb = bb_of_lambda(lambda);
bb_map = repmat(bb,1,1024);
pix = transpose(1:1024);

if showFig
    figure;
    plot(flipud(pix),spec)
    legend('spec pre R')
    title(sprintf('%s',file))
end

%inv_response_function = lflat./bb_map;
%inv_response_function = sum(inv_response_function,2);

inv_response_function = sum(lflat,2)./bb;

spec = spec./inv_response_function; % spec now in ph/sec/m2/sr/um

if showFig
    figure;
    plot(flipud(pix), 1./inv_response_function)
    title(sprintf('%s',file))
    legend('Response function')
    
    figure;
    plot(flipud(pix), spec)
    legend('spec processed')
    title(sprintf('%s',file))
    
    figure;
    plot(lambda, spec);
    xlabel('$\lambda$ [$\mu$m]', 'Interpreter', 'latex')
    ylabel('Intensity [ph/s/m$^2$/sr/$\mu$ m]')
end

%%
% ldata_plot = ldata;
% ldata_plot(ldata_plot > 150) = 150;
% dlmwrite(strcat(figpath, '\\data_lin_sat150.txt'), ldata_plot) % Save data for plotting later
%%

Mbti = dlmread('bandToIntervalAndNDF.txt');%Area of interest
isu = Mbti(Mbti(:,1) == trans,2);
ieu = Mbti(Mbti(:,1) == trans,3);

%% Overlap

mask = ones(size(spec)); %If you ever want to get rid of specific pixels

mask(lambda<isu)=0;
mask(lambda>ieu)=0;

Num = [0.4027,-2.0136,4.0273,-4.0273,2.0136,-0.4027];%filter
Den = [1.0000,-3.2160,4.3552,-3.0536,1.1003,-0.1622];

sspec=filtfilt(Num,Den,spec);
sspec(mask == 0) = 0;

if (filter=='H')
%% H filter processing    
    %% Noise cancelling 31 Q-branch
    pix = 1:1024;
    lambdaT = transpose(flipud(lambda));
    specT = transpose(flipud(spec)); % To align spec with pixel number again
    spec = flipud(spec); % to align spec with pixel number
    lambda = flipud(lambda); % Align lambda with pixel and spec
    %specT = transpose(spec);
    
    %q_31_peakrange = 110:133;
    %q_31_peakrange = 891:914;
    q_31_peakrange = 891:914;
    %Q_31_range = 100:160;
    %Q_31_range = 861:947;
    Q_31_range = 870:934;
    q31spec = spec; 
    fitq31spec = fit(transpose([pix(Q_31_range)]),q31spec(Q_31_range),'poly1');
    q31spec = q31spec - fitq31spec(pix);
    q31spec = q31spec + abs(min(q31spec(Q_31_range)));
    thresh = 0.05*mean(q31spec(Q_31_range));
    [~, locs, ~, p] = findpeaks(q31spec(q_31_peakrange), 'MinPeakWidth', 1,'MinPeakProminence',thresh);%,'MinPeakWidth', 1); %,'SortStr', 'descend');
    [~, loci] = max(p);
    q31i = locs(loci);
    %[~, q31i] = max(spec(q_31_peakrange));
    q31i = q31i + q_31_peakrange(1) - 1;
    %Q_31_spec = (q31i-5):(q31i+25);
    Q_31_spec = (q31i-25):(q31i+5);
    q31range = [q31i - 43, Q_31_spec(1) - 1 , Q_31_spec(end) + 1, q31i + 28];
    
%     % Find linear noise
%     if strcmp(file(1:6), 'NCyi28')
%         fit_type = 'poly1';
%     else
%         fit_type = 'poly3';
%     end
%     noise_per_pix_31 = fit(transpose([q31range(1):q31range(2) q31range(3):q31range(4)]),...
%         transpose([specT(q31range(1):q31range(2)) specT(q31range(3):q31range(4))]), fit_type);

    warning('off','all')
    noise_per_pix_31_3 = fit(transpose([q31range(1):q31range(2) q31range(3):q31range(4)]),...
        transpose([specT(q31range(1):q31range(2)) specT(q31range(3):q31range(4))]), 'poly3');
    
    noise_per_pix_31_2 = fit(transpose([q31range(1):q31range(2) q31range(3):q31range(4)]),...
        transpose([specT(q31range(1):q31range(2)) specT(q31range(3):q31range(4))]), 'poly2');
    warning('on','all')
%     nc_spec31 = spec - noise_per_pix_31(pix); % This is only valid in the 31 Q branch region!
%     nc_spec31h = transpose(nc_spec31);
    
    nc_spec31_2 = spec - noise_per_pix_31_2(pix); % This is only valid in the 31 Q branch region!
    nc_spec31h_2 = transpose(nc_spec31_2);
    
    nc_spec31_3 = spec - noise_per_pix_31_3(pix); % This is only valid in the 31 Q branch region!
    nc_spec31h_3 = transpose(nc_spec31_3);
    
    dlambda = abs(((lambda(Q_31_spec(end))-lambda(Q_31_spec(1)))/(length(Q_31_spec)-1)));
    counts31 = sum(nc_spec31_3(Q_31_spec))*dlambda;
    
    counts31_corr23 = mean([sum(nc_spec31_2(Q_31_spec))*dlambda, sum(nc_spec31_3(Q_31_spec))*dlambda]);
    
    stderror31 = std(transpose([nc_spec31h_3(q31range(1):q31range(2)) nc_spec31h_3(q31range(3):q31range(4))]))...
        /sqrt(length(transpose([Q_31_spec])));
    
    counts31_1 = counts31;
    counts31 = counts31_corr23; % Use mean of poly2 and poly3 fit to get the counts.
    
    spec_stderr31 = stderror31*length(Q_31_spec)*dlambda + (std([counts31,counts31_corr23])/sqrt(2));
    spec_relerr31 = spec_stderr31/counts31;
    
%     fprintf(' Total photons/sec/m2/sr from 31 = %d \n Relative error = %0.2f %% \n',counts31, spec_relerr31*100);
%  
    if showFig
        sprintf('Specific dark: %d, %s, 31: I_{3} = %g, I_{2+3} = %g',specificDark, S.PrimaryData.Keywords{eindex, 2},counts31_1,counts31)
        figure;
        plot(Q_31_range, specT(Q_31_range), Q_31_range, noise_per_pix_31_2(Q_31_range),Q_31_range, noise_per_pix_31_3(Q_31_range))
        title(sprintf('%s 31 Q-branch, rel.err = %0.2f %%',file,spec_relerr31*100))
        xlabel('Pix')
        ylabel('Intensity [a.u.]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([(q31range(2)+1) (q31range(2)+1)],y1, '--m', 'linewidth', 1)
        plot([(q31range(3)-1) (q31range(3)-1)],y1, '--m', 'linewidth', 1)
        hold off
        
        figure;
        plot(lambdaT(Q_31_range), specT(Q_31_range), lambda(Q_31_range), noise_per_pix_31_2(Q_31_range),lambdaT(Q_31_range), noise_per_pix_31_3(Q_31_range))
        title(sprintf('%s 31 Q-branch, rel.err = %0.2f %%',file,spec_relerr31*100))
        xlabel('$\lambda$ [$\mu$m ]')
        ylabel('Intensity [photons/s/m$^2$/sr/$\mu$m]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambdaT((q31range(2)+1)) lambdaT((q31range(2)+1))],y1, '--m', 'linewidth', 1)
        plot([lambdaT(q31range(3)-1) lambdaT(q31range(3)-1)],y1, '--m', 'linewidth', 1)
        hold off
    end

    %% Noise cancelling 42 Q-branch
    %Q_42_range = 350:439;
    Q_42_range = 586:682;
%     q_42_peakrange = 380:407;
%     thresh = 0.0005*mean(spec(Q_42_range));
%     q42spec = spec;
%     fitq42spec = fit(transpose([pix(Q_42_range)]),q42spec(Q_42_range),'poly1');
%     q42spec = q42spec - fitq42spec(pix);
%     [~, locs, ~, p] = findpeaks(q42spec(q_42_peakrange),'MinPeakProminence',thresh, 'MinPeakWidth', 0.6,'SortStr', 'descend'); % ,'MinPeakProminence',thresh);%,'MinPeakWidth', 1); %,'SortStr', 'descend');
%     [~, loci] = max(p);
%     q42i = locs(loci);
%     %[~, q31i] = max(spec(q_31_peakrange));
%     q42i = q42i + q_42_peakrange(1) - 1;
    q42i = q31i - 267;
    %q42i = q31i - 268;
    Q_42_spec = (q42i-25):(q42i+5);
    
    %q42range = [Q_42_range(1), Q_42_spec(1) - 1 , Q_42_spec(end) + 1, q42i + 42];
    %q42fitrange = [q42i-39, q42i-12, Q_42_spec(end) + 1, q42i + 42];
    
    q42range = [Q_42_range(1), Q_42_spec(1) - 1 , Q_42_spec(end) + 1, q42i + 40];
    q42fitrange = [q42i-44, Q_42_spec(1) - 1, q42i + 12, q42i + 40];

    lambdaT = transpose(lambda);
    specT = transpose(spec);
    
    warning('off','all')
    noise_per_pix_42 = fit(transpose([q42fitrange(1):q42fitrange(2) Q_42_spec(1) q42fitrange(3):q42fitrange(4)]),...
        transpose([specT(q42fitrange(1):q42fitrange(2)) specT(Q_42_spec(1)) specT(q42fitrange(3):q42fitrange(4))]), 'poly3');
    noise_per_pix_42_2 = fit(transpose([q42fitrange(1):q42fitrange(2) Q_42_spec(1) q42fitrange(3):q42fitrange(4)]),...
        transpose([specT(q42fitrange(1):q42fitrange(2)) specT(Q_42_spec(1)) specT(q42fitrange(3):q42fitrange(4))]), 'poly2');
    warning('on','all')
    
    nc_spec42 = spec - noise_per_pix_42(pix); % This is only valid in the 42 Q branch region!
    nc_spec42h = transpose(nc_spec42);
    nc_spec42_2 = spec - noise_per_pix_42_2(pix); % This is only valid in the 42 Q branch region!
    
    dlambda = abs(((lambda(Q_42_spec(end))-lambda(Q_42_spec(1)))/(length(Q_42_spec)-1)));
    counts42 = sum(nc_spec42(Q_42_spec))*dlambda;
    counts42_2 = sum(nc_spec42_2(Q_42_spec))*dlambda;
    
    counts42_23 = mean([counts42, counts42_2]); % Mean of counts using poly2 and poly3
    stderror42 = std(transpose([nc_spec42h(q42fitrange(1):q42fitrange(2)) nc_spec42h(q42fitrange(3):q42fitrange(4))]))...
       /sqrt(length(transpose([Q_42_spec])));
    
    spec_stderr42 = stderror42*length(Q_42_spec)*dlambda + abs(std([counts42_23, counts42])/sqrt(2));
    spec_relerr42 = spec_stderr42/counts42_23;
    
    if counts42_23 < 0 % Just sum the first peak, the rest should average to zero by white noise assumption
        counts42_23 = sum(nc_spec42(Q_42_spec(1):(Q_42_spec(1)+12)))*dlambda;
        spec_stderr42 = stderror42*length(Q_42_spec(1):(Q_42_spec(1)+12))*dlambda + abs(std([counts42_23, counts42])/sqrt(2));
        spec_relerr42 = spec_stderr42/counts42_23;
    end
   
    
    counts42_1 = counts42;
    counts42 = counts42_23;
    
    if showFig
        sprintf('Specific dark: %d, %s, 42: I_{3} = %g, I_{2+3} = %g',specificDark, S.PrimaryData.Keywords{eindex, 2},counts42_1,counts42_23)
        figure;
        plot(Q_42_range, specT(Q_42_range), Q_42_range, noise_per_pix_42(Q_42_range),Q_42_range, noise_per_pix_42_2(Q_42_range))
        title(sprintf('%s 42 Q-branch, rel.err = %0.2f %%',file,spec_relerr42*100))
        xlabel('Pix')
        ylabel('Intensity [a.u.]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([Q_42_spec(1) Q_42_spec(1)],y1, '--m', 'linewidth', 1)
        plot([Q_42_spec(end) Q_42_spec(end)],y1, '--m', 'linewidth', 1)
        hold off
        figure;
        plot(lambdaT(Q_42_range), specT(Q_42_range), lambdaT(Q_42_range), noise_per_pix_42(Q_42_range),lambdaT(Q_42_range), noise_per_pix_42_2(Q_42_range))
        title(sprintf('%s 42 Q-branch, rel.err = %0.2f %%',file,spec_relerr42*100))
        xlabel('$\lambda$ [$\mu$m ]')
        ylabel('Intensity [photons/s/m$^2$/sr/$\mu$m]') 
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambdaT(Q_42_spec(1)) lambdaT(Q_42_spec(1))],y1, '--m', 'linewidth', 1)
        plot([lambdaT(Q_42_spec(end)) lambdaT(Q_42_spec(end))],y1, '--m', 'linewidth', 1)
        hold off
    end

    %% Noise cancelling 53 Q-branch
    lambdaT = transpose(lambda);
    specT = transpose(spec);

    %Q_53_range = 610:710;
    Q_53_range = 317:380;
    q53i = q31i - 538;
    Q_53_spec = (q53i-25):(q53i+5);
    %q53fitrange = [q53i-23, q53i-4, q53i+24, q53i + 45];
    q53fitrange = [q53i-42, q53i-24, q53i+4, q53i + 20];
    %noise_per_lambda_53 = fit(transpose([lambdaT(618:651) lambdaT(673:698)]), transpose([specT(618:651) specT(673:698)]), 'poly1'); % This should be smarter
    
    warning('off','all')
    noise_per_pix_53 = fit(transpose([q53fitrange(1):q53fitrange(2) q53fitrange(3):q53fitrange(4)]),...
        transpose([specT(q53fitrange(1):q53fitrange(2)) specT(q53fitrange(3):q53fitrange(4))]), 'poly3');
    noise_per_pix_53_2 = fit(transpose([q53fitrange(1):q53fitrange(2) q53fitrange(3):q53fitrange(4)]),...
        transpose([specT(q53fitrange(1):q53fitrange(2)) specT(q53fitrange(3):q53fitrange(4))]), 'poly2');
    warning('on','all')
    nc_spec53 = spec - noise_per_pix_53(pix); % This is only valid in the 53 Q branch region!
    nc_spec53h = transpose(nc_spec53);
    
    nc_spec53_2 = spec - noise_per_pix_53_2(pix); % This is only valid in the 53 Q branch region!
    nc_spec53h_2 = transpose(nc_spec53_2);
    
    dlambda = abs(((lambda(Q_53_spec(end))-lambda(Q_53_spec(1)))/(length(Q_53_spec)-1)));
    counts53 = sum(nc_spec53(Q_53_spec))*dlambda;
    
    counts53_23 = mean([counts53, sum(nc_spec53_2(Q_53_spec))*dlambda]);
%     if counts42 < 0 % Just sum the first peak, the rest should average to zero by white noise assumption
%         counts42 = sum(nc_spec42(Q_42_spec(1):Q_42_spec(1)+12))*dlambda;
%     end
    
    stderror53 = std(transpose([nc_spec53h(q53fitrange(1):q53fitrange(2)) nc_spec53h(q53fitrange(3):q53fitrange(4))]))...
        /sqrt(length(transpose([Q_53_spec])));
    
    counts53_1 = counts53;
    counts53 = counts53_23;
    
    spec_stderr53 = stderror53*length(Q_53_spec)*dlambda + abs(std([counts53_23, counts53_1])/sqrt(2));
    spec_relerr53 = spec_stderr53/counts53;

    if showFig
        sprintf('Specific dark: %d, %s, 53: I_{3} = %g, I_{2+3} = %g',specificDark, S.PrimaryData.Keywords{eindex, 2},counts53_1,counts53_23)
        figure;
        plot(Q_53_range, specT(Q_53_range), Q_53_range, noise_per_pix_53(Q_53_range), Q_53_range, noise_per_pix_53_2(Q_53_range))
        title(sprintf('%s 53 Q-branch, rel.err = %0.2f %%',file,spec_relerr53*100))
        xlabel('Pix')
        ylabel('Intensity [a.u.]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([Q_53_spec(1) Q_53_spec(1)],y1, '--m', 'linewidth', 1)
        plot([Q_53_spec(end) Q_53_spec(end)],y1, '--m', 'linewidth', 1)
        hold off
        
        figure;
        plot(lambdaT(Q_53_range), specT(Q_53_range), lambdaT(Q_53_range), noise_per_pix_53(Q_53_range), lambdaT(Q_53_range), noise_per_pix_53_2(Q_53_range))
        title(sprintf('%s 53 Q-branch, rel.err = %0.2f %%',file,spec_relerr53*100))
        xlabel('$\lambda$ [$\mu$m ]')
        ylabel('Intensity [photons/s/m$^2$/sr/$\mu$m]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambdaT(Q_53_spec(1)) lambdaT(Q_53_spec(1))],y1, '--m', 'linewidth', 1)
        plot([lambdaT(Q_53_spec(end)) lambdaT(Q_53_spec(end))],y1, '--m', 'linewidth', 1)
        hold off       
    end
    
    %% Noise cancelling 64 Q-branch 
    lambdaT = transpose(lambda);
    specT = transpose(spec);
    %Q_64_range = 920:990;
    Q_64_range = 43:104;
    q64i = q31i - 821;
    %Q_64_spec = (q64i-6):(q64i+26);
    Q_64_spec = (q64i-26):(q64i+6);
    %q64fitrange = [q64i-29, q64i-6, q64i+20, q64i + 42];
    q64fitrange = [q64i-47, q64i-25, q64i+5, q64i + 33];
    
%     if strcmp(file(1:6), 'NCsa10')
%         fit_type = 'poly1';
%     else
%         fit_type = 'poly3';
%     end
    warning('off','all')
    noise_per_pix_64 = fit(transpose([q64fitrange(1):q64fitrange(2) q64fitrange(3):q64fitrange(4)]),...
        transpose([specT(q64fitrange(1):q64fitrange(2)) specT(q64fitrange(3):q64fitrange(4))]), 'poly3');
    noise_per_pix_corr2 = fit(transpose([q64fitrange(1):q64fitrange(2) q64fitrange(3):q64fitrange(4)]),...
        transpose([specT(q64fitrange(1):q64fitrange(2)) specT(q64fitrange(3):q64fitrange(4))]), 'poly2');
    noise_per_pix_corr1 = fit(transpose([q64fitrange(1):q64fitrange(2) q64fitrange(3):q64fitrange(4)]),...
        transpose([specT(q64fitrange(1):q64fitrange(2)) specT(q64fitrange(3):q64fitrange(4))]), 'poly1');
    warning('on','all')
    nc_spec64 = spec - noise_per_pix_64(pix); % This is only valid in the 53 Q branch region!
    nc_spec64h = transpose(nc_spec64);
    
    nc_spec64_corr2 = spec - noise_per_pix_corr2(pix);
    nc_spec64_corr1 = spec - noise_per_pix_corr1(pix);
    
    noise_fit_corr = fit(transpose([Q_64_spec(1) Q_64_spec(end)]), transpose([specT(Q_64_spec(1)) specT(Q_64_spec(end))]),'poly1');
    nc_spec64_corr = spec - noise_fit_corr(pix);
    
    dlambda = abs(((lambda(Q_64_spec(end))-lambda(Q_64_spec(1)))/(length(Q_64_spec)-1)));
    counts64 = sum(nc_spec64(Q_64_spec))*dlambda; % With poly3
    counts64_corr = mean([counts64, sum(nc_spec64_corr(Q_64_spec))*dlambda]);
    counts64_corr2 = mean([counts64, sum(nc_spec64_corr(Q_64_spec))*dlambda,sum(nc_spec64_corr2(Q_64_spec))*dlambda]);
    
    counts64_corr3 = mean([counts64,sum(nc_spec64_corr2(Q_64_spec))*dlambda]); % Poly2 + poly3
    
    counts64_fit123 = dlambda*mean([sum(nc_spec64(Q_64_spec)), sum(nc_spec64_corr1(Q_64_spec)), sum(nc_spec64_corr2(Q_64_spec))]);
    
    stderror64 = std(transpose([nc_spec64h(q64fitrange(1):q64fitrange(2)) nc_spec64h(q64fitrange(3):q64fitrange(4))]))...
        /sqrt(length(transpose([Q_64_spec])));
    
    spec_stderr64 = stderror64*length(Q_64_spec)*dlambda + abs(std([counts64_corr3, counts64])/sqrt(2));
    spec_relerr64 = spec_stderr64/counts64_corr3;
%     if str2double(expMode{3}) < 3 % If number of ramps < 3 the noise is too heavy for proper 64
%         counts64 = sum(nc_spec64((q64i-5):(q64i+8)))*dlambda; % Sum just the two first peaks
%         counts64_corr = mean([counts64, sum(nc_spec64_corr((q64i-5):(q64i+8)))*dlambda]);
%         counts64_corr2 = mean([counts64, sum(nc_spec64_corr((q64i-5):(q64i+8)))*dlambda,sum(nc_spec64_corr2((q64i-5):(q64i+8)))*dlambda]);
%         counts64_corr3 = mean([counts64, sum(nc_spec64_corr2((q64i-5):(q64i+8)))*dlambda]);
%     end
    if counts64_fit123 < 0 % Just sum the first peak, the rest should average to zero by white noise assumption
        %counts64 = sum(nc_spec64((q64i-5):(q64i+5)))*dlambda;
        counts64_corr = mean([counts64, sum(nc_spec64_corr((q64i-5):(q64i+5)))*dlambda]);
        counts64_corr2 = mean([counts64, sum(nc_spec64_corr((q64i-5):(q64i+5)))*dlambda,sum(nc_spec64_corr2((q64i-5):(q64i+8)))*dlambda]);
        counts64_corr3 = mean([counts64, sum(nc_spec64_corr2((q64i-5):(q64i+5)))*dlambda]);
        counts64_fit123 = dlambda*mean([sum(nc_spec64((q64i-5):(q64i+5))),...
            sum(nc_spec64_corr1((q64i-5):(q64i+5))), sum(nc_spec64_corr2((q64i-5):(q64i+5)))]);
    end
    
    if counts64_corr3 < 0
        counts64_corr3 = mean([sum(nc_spec64((q64i-8):(q64i+5)))*dlambda,sum(nc_spec64_corr2((q64i-8):(q64i+5)))*dlambda]);
        spec_stderr64 = stderror64*length((q64i-8):(q64i+5))*dlambda + abs(std([counts64_corr3, counts64])/sqrt(2));
        spec_relerr64 = spec_stderr64/counts64_corr3;
    end
    if counts64_corr3 < 0
        counts64_corr3 = mean([sum(nc_spec64((q64i-5):(q64i+5)))*dlambda,sum(nc_spec64_corr2((q64i-5):(q64i+5)))*dlambda]);
        spec_stderr64 = stderror64*length((q64i-5):(q64i+5))*dlambda + abs(std([counts64_corr3, counts64])/sqrt(2));
        spec_relerr64 = spec_stderr64/counts64_corr3;
    end
    
    counts64_1 = counts64;
    counts64 = counts64_corr3;

    if showFig
        sprintf('Specific dark: %d, %s, 64: I_{uncorrected} = %g, I_{2+3} = %g',specificDark, S.PrimaryData.Keywords{eindex, 2},counts64_1,counts64_corr3)
        figure;
        plot(Q_64_range, specT(Q_64_range), Q_64_range, noise_per_pix_64(Q_64_range),...
            Q_64_range, noise_per_pix_corr2(Q_64_range))
        title(sprintf('%s %s 64 Q-branch, rel.err = %0.2f %%',file,S.PrimaryData.Keywords{eindex, 2},spec_relerr64*100))
        xlabel('Pix')
        ylabel('Intensity [a.u.]')  
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([Q_64_spec(1) Q_64_spec(1)],y1, '--m', 'linewidth', 1)
        plot([Q_64_spec(end) Q_64_spec(end)],y1, '--m', 'linewidth', 1)
        hold off
        
        figure;
        plot(lambdaT(Q_64_range), specT(Q_64_range), lambdaT(Q_64_range), noise_per_pix_64(Q_64_range),...
            lambdaT(Q_64_range), noise_per_pix_corr2(Q_64_range))
        title(sprintf('%s %s 64 Q-branch, rel.err = %0.2f %%',file,S.PrimaryData.Keywords{eindex, 2},spec_relerr64*100))
        xlabel('$\lambda$ [$\mu$m]')
        ylabel('Intensity [photons/s/m$^2$/sr/$\mu$m]')
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambdaT(Q_64_spec(1)) lambdaT(Q_64_spec(1))],y1, '--m', 'linewidth', 1)
        plot([lambdaT(Q_64_spec(end)) lambdaT(Q_64_spec(end))],y1, '--m', 'linewidth', 1)
        hold off
    end
    
%     if exposureTime < 20
%         counts64 = 0;
%     end
    
    %% Save data
    spec = flipud(spec); % Flip spec back to wavelength alignment
    sM = [lambda, spec, sigma, sspec, sflat, wavenumber];%, nc_spec31, nc_spec42, nc_spec53, nc_spec64];
    dlmwrite(strcat(figpath, sprintf('%s-spectrum.txt', filter)), sM) % Birk's pref way of saving

elseif filter == 'K'
    %% Noise cancelling 97 Q-branch
    Q_97_range = 536:588;
    Q_97_spec = 545:569;
    q97range = [Q_97_range(1), Q_97_spec(1) - 1 , Q_97_spec(end) + 1, Q_97_range(end)];
    
    lambdaT = transpose(lambda);
    specT = transpose(spec);
    noise_per_lambda_97 = fit(transpose([lambdaT(q97range(1):q97range(2)) lambdaT(q97range(3):q97range(4))]),...
                                         transpose([specT(q97range(1):q97range(2)) specT(q97range(3):q97range(4))]), 'poly2');
    figure;
    plot(lambda(Q_97_range), spec(Q_97_range), lambda(Q_97_range), noise_per_lambda_97(lambda(Q_97_range)))
    title('noise fit for 97 Q-branch')
    xlabel('\lambda [{\mu}m]')
    ylabel('Intensity [ADU]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([lambda(q97range(2)+1) lambda(q97range(2)+1)],y1, '--m', 'linewidth', 1)
    plot([lambda(q97range(3)-1) lambda(q97range(3)-1)],y1, '--m', 'linewidth', 1)
    hold off
    
else
%% Other filters
    %% Save processed data
    sM = [lambda, spec, sigma, sspec, sflat, wavenumber];
    dlmwrite(strcat(figpath, sprintf('%s-spectrum.txt', filter)), sM) % Birk's pref way of saving
    % save(strcat(path,sprintf('%d-',trans),file,'.txt'),'sM','-ascii');
end 

end