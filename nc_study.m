%% Read in data

wd = 'C:\Users\birke\Documents\prosjekt\ReduceSingleNotfileBirk\'; % working directory
%wd = 'D:\NOTArchive\';
filter = 'H';

file = ['NCwb191105'; 'NCwb191106'; 'NCwb191107'; 'NCwb191108']; % H
%file = ['NCwb191105'];
%file = ['NCwb191101'; 'NCwb191102'; 'NCwb191103'; 'NCwb191104']; % K
%file = ['NCwb191101'];
%file = ['NCwb191109'; 'NCwb191110'; 'NCwb191111'; 'NCwb191112']; % J

%file = ['NCza230395';'NCza230395'; 'NCza230396'; 'NCza230397';'NCza230398';'NCza230399';'NCza230400';'NCza230401'];
    
% day = 'NCyi28';
% file_cell = findHfiles(day);
% [r,c] = size(file_cell);
% file = [];
% for i = 1:r
%     file = [file; file_cell{i}];
% end

spec_struct = struct([]);
for i=1:size(file)
    SinglePrepare_dev(file(i,:),'none','none',false);
    figpath(i, :) = [wd, 'data\', sprintf('%s-%s\\',filter, file(i,:))]; % To save/read figures for debug/study
    M = dlmread(strcat(figpath(i,:),sprintf('%s-spectrum.txt', filter)));
    spec_struct(1).spectrum(:,i) = [M(:,2)];
end

%M = dlmread(strcat(figpath,sprintf('%s-spectrum.txt', filter)));
lambda = M(:,1);
spec = M(:,2);
%nc_spec42 = M(:,3);
sigma = M(:,3);
sspec = M(:,4);
sflat = M(:,5);
wavenumber = M(:,6);

%figure;
%plot(spec_struct.spectrum(:,2))

added_spec = sum(spec_struct.spectrum(:,:), 2); % Adds all the spectras together

% figure;
% plot(wavenumber, added_spec)

if filter == 'K'
    %% branch
end
if filter == 'K'
    %% 97 Q branch
    Q_97_range = 536:588;
    Q_97_spec = 545:569;
    %q97range = [536, 545 - 1 , 569 + 1, 588];
    q97range = [Q_97_range(1), Q_97_spec(1) - 1 , Q_97_spec(end) + 1, Q_97_range(end)];
    
    lambdaT = transpose(lambda);
    specT = transpose(added_spec);
    noise_per_lambda_97 = fit(transpose([lambdaT(q97range(1):q97range(2)) lambdaT(q97range(3):q97range(4))]),...
                                         transpose([specT(q97range(1):q97range(2)) specT(q97range(3):q97range(4))]), 'poly1');
                                     
    nc_added_spec97 = added_spec - noise_per_lambda_97(lambda); % This is only valid in the 97 Q branch region!
    
    stderror97 = std(transpose([specT(q97range(1):q97range(2)) specT(q97range(3):q97range(4))]))...
        /sqrt(length(transpose([specT(Q_97_spec)])));
    
    error_estimate97 = stderror97*length(Q_97_spec);
    
    counts = sum(nc_added_spec97((Q_97_spec)));
    
    figure;
    plot(lambda(Q_97_range), added_spec(Q_97_range), lambda(Q_97_range), noise_per_lambda_97(lambda(Q_97_range)))
    title('noise fit of added spectra for 97 Q-branch')
    xlabel('\lambda [{\mu}m]')
    ylabel('Intensity [ADU]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([lambda(q97range(2)+1) lambda(q97range(2)+1)],y1, '--m', 'linewidth', 1)
    plot([lambda(q97range(3)-1) lambda(q97range(3)-1)],y1, '--m', 'linewidth', 1)
    hold off

end

if filter == 'H'
    %% 31 Q branch
%     Q1_31 = 0;
%     Q_31_range = 105:161;
%     Q_31_spec = 114:133;
%     q31range = [105, 114, 133, 161];
    
    added_spec = flipud(added_spec);
    pix = 1:1024;
    q_31_peakrange = 110:140;
    q_31_peakrange = 891:914;
    Q_31_range = 100:160;
    Q_31_range = 861:934;
    [~, q31i] = max(added_spec(q_31_peakrange));
    q31i = q31i + q_31_peakrange(1) - 1;
    Q_31_spec = (q31i-25):(q31i+5);
    q31range = [Q_31_range(1), Q_31_spec(1) - 1 , Q_31_spec(end) + 1, Q_31_range(end)];

    lambdaT = transpose(lambda);
    specT = transpose(added_spec);
    %noise_per_lambda_31 = fit(transpose([lambdaT(105:113) lambdaT(134:161)]), transpose([specT(105:113) specT(134:161)]), 'poly1');
    
    noise_per_pix_31 = fit(transpose([q31range(1):q31range(2) q31range(3):q31range(4)]),...
        transpose([specT(q31range(1):q31range(2)) specT(q31range(3):q31range(4))]), 'poly3');
    
    nc_added_spec31 = added_spec - noise_per_pix_31(pix); % This is only valid in the 31 Q branch region!
    
    stderror31 = std(transpose([specT(q31range(1):q31range(2)) specT(q31range(3):q31range(4))]))...
        /sqrt(length(transpose([specT(Q_31_spec)])));

    figure;
    hold on;
    plot(Q_31_range, added_spec(Q_31_range)) 
    plot(Q_31_range, noise_per_pix_31(Q_31_range))
    title('noise fit of added spectra for 31 Q-branch')
    xlabel('Pix]')
    ylabel('Intensity [a.u.]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([(q31range(2)+1) (q31range(2)+1)],y1, '--m', 'linewidth', 1)
    plot([(q31range(3)-1) (q31range(3)-1)],y1, '--m', 'linewidth', 1)
    hold off


    %% 42 Q branch

    Q_42_range = 586:667;
%     q_42_peakrange = 380:410;
%     thresh = 0.01*mean(spec(Q_42_range));
%     [~, locs, ~, p] = findpeaks(spec(q_42_peakrange),'MinPeakProminence',thresh, 'MinPeakWidth', 1,'SortStr', 'descend'); % ,'MinPeakProminence',thresh);%,'MinPeakWidth', 1); %,'SortStr', 'descend');
%     %[~, loci] = max(p);
%     q42i = locs(1);
%     %[~, q31i] = max(spec(q_31_peakrange));
%     q42i = q42i + q_42_peakrange(1) - 1;
    q42i = q31i - 267;
    Q_42_spec = (q42i-25):(q42i+5);
    
    q42range = [Q_42_range(1), Q_42_spec(1) - 1 , Q_42_spec(end) + 1, q42i + 42];
    q42fitrange = [q42i-39, q42i-12, Q_42_spec(end) + 1, q42i + 42];
    q42fitrange = [q42i-42, Q_42_spec(1) - 1, q42i + 12, q42i + 39];

    lambdaT = transpose(lambda);
    specT = transpose(added_spec);
    
    noise_per_pix_42 = fit(transpose([q42fitrange(1):q42fitrange(2) q42fitrange(3):q42fitrange(4)]),...
        transpose([specT(q42fitrange(1):q42fitrange(2)) specT(q42fitrange(3):q42fitrange(4))]), 'poly3');

    nc_added_spec42 = added_spec - noise_per_pix_42(pix); % This is only valid in the 42 Q branch region!

    figure;
    hold on;
    plot(Q_42_range, added_spec(Q_42_range)) 
    plot(Q_42_range, noise_per_pix_42(Q_42_range))
    title('noise fit of added spectra for 42 Q-branch')
    xlabel('Pix]')
    ylabel('Intensity [a.u.]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([Q_42_spec(1) Q_42_spec(1)],y1, '--m', 'linewidth', 1)
    plot([Q_42_spec(end) Q_42_spec(end)],y1, '--m', 'linewidth', 1)
    hold off

    %% 53 Q branch

    Q_53_range = 317:380;
    q53i = q31i - 538;
    Q_53_spec = (q53i-25):(q53i+5);
    %q53fitrange = [q53i-20, q53i-4, q53i+24, q53i + 42];
    q53fitrange = [q53i-42, q53i-24, q53i+4, q53i + 20];

    lambdaT = transpose(lambda);
    specT = transpose(added_spec);
    %noise_per_lambda_53 = fit(transpose([lambdaT(618:651) lambdaT(673:698)]), transpose([specT(618:651) specT(673:698)]), 'poly1'); % This should be smarter
    noise_per_pix_53 = fit(transpose([q53fitrange(1):q53fitrange(2) q53fitrange(3):q53fitrange(4)]),...
        transpose([specT(q53fitrange(1):q53fitrange(2)) specT(q53fitrange(3):q53fitrange(4))]), 'poly3');
    nc_added_spec53 = added_spec - noise_per_pix_53(pix); % This is only valid in the 53 Q branch region!


    figure;
    plot(Q_53_range, added_spec(Q_53_range), Q_53_range, noise_per_pix_53(Q_53_range))
    title('noise fit of added spectra for 53 Q-brach')
    xlabel('pix')
    ylabel('Intensity [ADU]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([Q_53_spec(1) Q_53_spec(1)],y1, '--m', 'linewidth', 1)
    plot([Q_53_spec(end) Q_53_spec(end)],y1, '--m', 'linewidth', 1)
    hold off
    
    %% 64 branch
    
    Q_64_range = 43:104;
    q64i = q31i - 821;
    Q_64_spec = (q64i-26):(q64i+5);
    %q64fitrange = [q64i-29, q64i-4, q64i+26, q64i + 42];
    q64fitrange = [q64i-42, q64i-26, q64i+4, q64i + 29];

    lambdaT = transpose(lambda);
    specT = transpose(added_spec);
    %noise_per_lambda_53 = fit(transpose([lambdaT(618:651) lambdaT(673:698)]), transpose([specT(618:651) specT(673:698)]), 'poly1'); % This should be smarter
    noise_per_pix_64 = fit(transpose([q64fitrange(1):q64fitrange(2) q64fitrange(3):q64fitrange(4)]),...
        transpose([specT(q64fitrange(1):q64fitrange(2)) specT(q64fitrange(3):q64fitrange(4))]), 'poly3');
    nc_added_spec64 = added_spec - noise_per_pix_64(pix); % This is only valid in the 53 Q branch region!


    figure;
    plot(Q_64_range, added_spec(Q_64_range), Q_64_range, noise_per_pix_64(Q_64_range))
    title('noise fit of added spectra for 64 Q-brach')
    xlabel('pix')
    ylabel('Intensity [ADU]')

    y1=get(gca,'ylim');
    hold on % Vertical lines to visualize integration area
    plot([Q_64_spec(1) Q_64_spec(1)],y1, '--m', 'linewidth', 1)
    plot([Q_64_spec(end) Q_64_spec(end)],y1, '--m', 'linewidth', 1)
    hold off
end


