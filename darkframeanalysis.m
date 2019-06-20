%% Dark frame analysis. How do they very with time?
wd = 'D:\NOTArchive\'; % working directory
master = importtxt('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NOTCamOverviewMaster.txt',',',2);
%%
RUN = true;
if RUN
%darkframes = cell(1,1231);
darkframes = cell(1,1231);
dexptimes = cell(1,1231);
ddate = cell(1,1231);
dspecs_unproc = cell(1,1231);
dspecs_proc = cell(1,1231);
%dark_data_raw = struct('Data_raw',[], 'Data_proc', []);
dark_data_raw = [zeros(1024,1024);[]];
dark_data_proc = [];
d = 1;
f = waitbar(0, 'Working!!');
for line = 1:size(master)
    expmode = split(master{line,4});
    if strcmp(expmode{1}, 'dframes')
            dfile = master{line,1}(1:10);
            dexptimes{d} = str2double(expmode{2})*str2double(expmode{3});
            darkframes{d} = dfile;
            Y = master{line,1}(3) + 1894;
            M = master{line,1}(4) - 96;
            D = str2double(master{line,1}(5:6));
            ddate{d} = datetime(Y, M, D);
            dpath = [wd, dfile(1:6), '\\'];
            dark = fitsread(strcat(dpath,dfile,'.fits'),'Image');
            dspecs_unproc{d} = sum(dark,2)/dexptimes{d};
            %dark_data_raw(:,:,d) = dark/dexptimes{d};
            bad = fitsread('bad_zero_sci.fits');
            bad(bad > 0)=1;
            dark(bad == 1) = 10001;
            %3 smooth over 10x10 pixels (median) -> dark
            hot = [handleBackground(dark(1:512,1:512)),handleBackground(dark(1:512,513:1024));handleBackground(dark(513:1024,1:512)),handleBackground(dark(513:1024,513:1024))];
            hot = [handleBackground2(hot(1:512,1:512)),handleBackground2(hot(1:512,513:1024));handleBackground2(hot(513:1024,1:512)),handleBackground2(hot(513:1024,513:1024))];
            dark(:,512)=0;
            dark(bad == 1) = 0;
            dark = dark - hot;
            %%
            hot(hot>550*dexptimes{d}/10.8)=0;
            %dark(dark<-10000)=0;
            hotaverage = mean(mean(hot));
            hotsigma = std2(hot);
            hot(abs(hot)>abs(hotaverage)+6*hotsigma)=-100001;%0.06
            hot(hot>=-100000)=1;
            hot(hot<-100000)=0;
            %disp(sum(sum(abs(hot),1))/1024/1024);
            dark(hot == 0) = 0;
            dspecs_proc{d} = sum(dark,2)/dexptimes{d};
            %dark_data_proc(:,:,d) = dark/dexptimes{d};
            waitbar(d/1231, f, 'Scanning darkframes...');
%             if d == 410
%                 masterDark_raw410 = mean(dark_data.Data_raw(:,:,:),3);
%                 %masterDark_proc615 = mean(dark_data.Data_proc(:,:,:),3);
%             elseif d == 820
            d = d + 1;
    end
end
close(f)
end

runMasterDark = true;
dexptimes = cell(1,1231);
dark_data_raw = [zeros(1024,1024, 412)];
dark_data_proc = [];
d = 1;
k = 1;
if runMasterDark
f = waitbar(0, 'Working!!');
for line = 1:size(master)
    expmode = split(master{line,4});
    if strcmp(expmode{1}, 'dframes')
            dfile = master{line,1}(1:10);
            dexptimes{d} = str2double(expmode{2})*str2double(expmode{3});
            dpath = [wd, dfile(1:6), '\\'];
            dark = fitsread(strcat(dpath,dfile,'.fits'),'Image');
            dark_data_raw(:,:,d) = dark/dexptimes{d};
            waitbar(k/1231, f, 'Preparing master darks');
            if k == 409
                masterDark_raw409 = median(dark_data_raw(:,:,1:409),3);
                d = 0;
            elseif k == 819
                masterDark_raw410 = median(dark_data_raw(:,:,1:410),3);
                d = 0;
            elseif k == 1231
                masterDark_raw412 = median(dark_data_raw(:,:,1:412),3);
                d = 0;
            end
            d = d + 1;
            k = k + 1;
    end
end
clear dark_data_raw;
close(f)
end

darkframes_copy = darkframes;
dexptimes_copy = dexptimes;
ddate_copy = ddate;
dspecs_unproc_copy = dspecs_unproc;
dspecs_proc_copy = dspecs_proc;

darks_wb19 = [dspecs_proc{495:498}];

darks_ze21 = [dspecs_proc{1160:1190}];
darks_ze22 = [dspecs_proc{1198:1211}];
darks_zf20 = [dspecs_proc{1216:1220}];
darks_zf21 = [dspecs_proc{1224:1231}];

darks_wb19u = [dspecs_unproc{495:498}];

darks_ze21u = [dspecs_unproc{1160:1190}];
darks_ze22u = [dspecs_unproc{1198:1211}];
darks_zf20u = [dspecs_unproc{1216:1220}];
darks_zf21u = [dspecs_unproc{1224:1231}];

masterDark_spec_raw = [mean([dspecs_unproc{:}],2)];
masterDark_raw = [];

masterDark_raw(:,:,1) = masterDark_raw409*409/412;
masterDark_raw(:,:,2) = masterDark_raw410*410/410;
masterDark_raw(:,:,3) = masterDark_raw412;

masterDark_raw = median(masterDark_raw,3);
%dlmwrite('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\masterDarks\masterDarks_raw_median.txt', masterDark_raw)
masterDark_mean = dlmread('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\masterDarks\masterDark_raw.txt');

figure;
plot([sum(masterDark_raw,2), sum(masterDark_mean,2)])
legend('median', 'mean')

figure;
pix = 1:1024;
plot(pix, [median(darks_wb19u,2), median(darks_ze21u,2), median(darks_ze22u,2), median(darks_zf20u,2), median(darks_zf21u,2)], 'b')
hold on
plot(pix, masterDark_spec_raw, '-r')

dark_intensity = zeros(1,1231);
for i = 1:1231
    intensity = sum([dspecs_proc{i}]);
    dark_intensity(i) = intensity;
end

figure;
plot([ddate{:}], dark_intensity, '.-')
legend('Darkframe intensity in sum of counts/s per pixel')









% for i = 1:length(darkframes)-2
%     if i == 1
%         darkframes_copy{1} = 9999;
%         dexptimes_copy{1} = 9999;
%         ddate_copy{1} = 9999;
%         dspecs_unproc_copy{1} = 9999;
%         dspecs_proc_copy{1} = 9999;
%     else
%         currentNumber = str2double(darkframes{i}(7:10));
%         nextNumber = str2double(darkframes{i+1}(7:10));
%         nextNextNumber = str2double(darkframes{i+2}(7:10));
%         if (nextNumber ~= (currentNumber + 1)) && ((nextNumber+1) == (nextNextNumber))
%             darkframes_copy{i+1} = 9999;
%             dexptimes_copy{i+1} = 9999;
%             ddate_copy{i+1} = 9999;
%             dspecs_unproc_copy{i+1} = 9999;
%             dspecs_proc_copy{i+1} = 9999;
%         end
%     end
% end


% darkframes_copy2 = {};
% dexptimes_copy2 = {};
% ddate_copy2 = {};
% dspecs_unproc_copy2 = {};
% dspecs_proc_copy2 = {};
% 
% new_index = 1;
% for i = 1:length(darkframes_copy)
%     if darkframes_copy{i} ~= 9999
%         darkframes_copy2{new_index} = darkframes_copy{i};
%         dexptimes_copy2{new_index} = dexptimes_copy{i};
%         ddate_copy2{new_index} = ddate_copy{i};
%         dspecs_unproc_copy2{new_index} = dspecs_unproc_copy{i};
%         dspecs_proc_copy2{new_index} = dspecs_proc_copy{i};
%         new_index = new_index + 1;
%     end
% end
% 
%     
% 
