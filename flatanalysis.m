%% Flatfields analysis. How do they vary with time?
wd = 'D:\NOTArchive\'; % working directory
master = importtxt('C:\Users\birke\Documents\master\scripting\NOTCamOverviewMaster.txt',',',2);
set(0,'defaulttextinterpreter','latex')
filter = 'J';
%%
if true
dfile = 'NCzc250213';
dpath = [wd, dfile(1:6), '\\'];
dark = fitsread(strcat(dpath,dfile,'.fits'),'Image');
bad = fitsread('bad_zero_sci.fits');
%3 smooth over 10x10 pixels (median) -> dark
hot = [handleBackground(dark(1:512,1:512)),handleBackground(dark(1:512,513:1024));handleBackground(dark(513:1024,1:512)),handleBackground(dark(513:1024,513:1024))];
hot = [handleBackground2(hot(1:512,1:512)),handleBackground2(hot(1:512,513:1024));handleBackground2(hot(513:1024,1:512)),handleBackground2(hot(513:1024,513:1024))];
dark(:,512)=0;
dark(bad == 1) = 0;
%dark = dark - hot;
hot(hot>550*600/10.8)=0;
%dark(dark<-10000)=0;
hotaverage = mean(mean(hot));
hotsigma = std2(hot);
hot(abs(hot)>abs(hotaverage)+6*hotsigma)=-100001;%0.06
hot(hot>=-100000)=1;
hot(hot<-100000)=0;
%disp(sum(sum(abs(hot),1))/1024/1024);
dark(hot == 0) = -1;
dark = dark/600;
%darkframes = cell(1,1231);
% darkframes = cell(1,1231);
% dexptimes = cell(1,1231);
% ddate = cell(1,1231);
% dspecs_unproc = cell(1,1231);
% dspecs_proc = cell(1,1231);
flatfields = {};
fexptimes = {};
fdate = {};
fspecs_unproc = {};
fspecs_proc = {};
fslit = [];
flat_data = struct('Data_raw',[], 'Data_proc', []);
f = 1;
%w = waitbar(0, 'Working!!');
for line = 1:size(master)
    expmode = split(master{line,4});
    if contains(master{line,3}, 'Halogen') && ~strcmp(master{line,3}, 'Halogen 1') && strcmp(master{line,2},filter)
        slit = 1;
        ffile = master{line,1}(1:10);
        flatfields{f} = ffile;
        fexptime = split(master{line, 4}); % {'expose', '7.0'}
        fexptime = str2double(fexptime(2));
        fexptimes{f} = fexptime;
        Y = master{line,1}(3) + 1894;
        M = master{line,1}(4) - 96;
        D = str2double(master{line,1}(5:6));
        fdate{f} = datetime(Y, M, D);
        fpath = [wd, ffile(1:6), '\\'];
        flat = fitsread(strcat(fpath,ffile,'.fits'),'Image');
        fspecs_unproc{f} = sum(flat,2)/fexptime;
        flat = flat/fexptime;
        flat_data.Data_raw(:,:,f) = flat;
        %substract zero -> put to high number
        flat = flat - dark;
        flat(bad == 1) = 0;%1000001;
        %substract hot -> put to high number
        flat(hot == 0) = 0;%1000001;
        %normalize that the average pixel is 1 in the middle (corners are gonna explode)
        %flat = flat/mean(mean(flat)); % Lets try not doing this
        flat(flat==0)=1;
        flat_data.Data_proc(:,:,f) = flat;
        fspecs_proc{f} = sum(flat,2);
        [~,S] = getFitsData(ffile);
        k = S.PrimaryData.Keywords;
        [aindex1, ~] = find(strcmp(k,'APERTUR'));
        [aindex2, ~] = find(strcmp(k,'NCAPRNM'));
        if ~isempty(aindex1)
            if ~strcmp(k{aindex1,2}, '128mu WF slit')
                slit = 0;
            end
        elseif ~isempty(aindex2)
            if ~strcmp(k{aindex2,2}, '128mu WF slit')
                slit = 0;
            end
        end
        fslit(f) = slit;
        f = f + 1;
    end
end
%close(w)
end

%% Clean flats 
if strcmp(filter, 'H')
    bad_flats = [3,6,21,24,28,29,30,51];
    flat_data.Data_raw(:,:,[bad_flats]) = [];
    flat_data.Data_proc(:,:,[bad_flats]) = [];
    fspecs_proc([bad_flats]) = []; 
    fspecs_unproc([bad_flats]) = []; 
    fdate([bad_flats]) = [];
    fslit([bad_flats]) = [];
elseif strcmp(filter, 'K')
    bad_flats = [1,19,20,45,46,49,50,51,54,55,68,69,70,73,78,79,92,93,94,95,96,97,98,99,104,105,114,115,120,121];
    flat_data.Data_raw(:,:,[bad_flats]) = [];
    flat_data.Data_proc(:,:,[bad_flats]) = [];
    fspecs_proc([bad_flats]) = []; 
    fspecs_unproc([bad_flats]) = []; 
    fdate([bad_flats]) = [];
    fslit([bad_flats]) = [];
elseif strcmp(filter, 'J')
    
end

masterFlat_raw = median([flat_data.Data_raw(:,:,:)],3);
dlmwrite('C:\Users\birke\Documents\master\scripting\masterFlats\J_masterFlat_raw_median.txt', masterFlat_raw)
% mean2009 = mean([fspecs_proc{1}],2);
% mean2010 = mean([fspecs_proc{2}],2);
% mean2011 = mean([fspecs_proc{3:9}],2);
% mean2012 = mean([fspecs_proc{10:30}],2);
% mean2013 = mean([fspecs_proc{31:38}],2);
% mean2014 = mean([fspecs_proc{39:46}],2);
% mean2015 = mean([fspecs_proc{47:73}],2);
% mean2016 = mean([fspecs_proc{74:79}],2);
% 
% mean2009u = mean([fspecs_unproc{1}],2);
% mean2010u = mean([fspecs_unproc{2}],2);
% mean2011u = mean([fspecs_unproc{3:9}],2);
% mean2012u = mean([fspecs_unproc{10:30}],2);
% mean2013u = mean([fspecs_unproc{31:38}],2);
% mean2014u = mean([fspecs_unproc{39:46}],2);
% mean2015u = mean([fspecs_unproc{47:73}],2);
% mean2016u = mean([fspecs_unproc{74:79}],2);

% mean2009 = median([fspecs_proc{1}],2);
% mean2010 = median([fspecs_proc{2}],2);
% mean2011 = median([fspecs_proc{3:9}],2);
% mean2012 = median([fspecs_proc{10:30}],2);
% mean2013 = median([fspecs_proc{31:38}],2);
% mean2014 = median([fspecs_proc{39:46}],2);
% mean2015 = median([fspecs_proc{47:73}],2);
% mean2016 = median([fspecs_proc{74:79}],2);
% 
% mean2009u = median([fspecs_unproc{1}],2);
% mean2010u = median([fspecs_unproc{2}],2);
% mean2011u = median([fspecs_unproc{3:9}],2);
% mean2012u = median([fspecs_unproc{10:30}],2);
% mean2013u = median([fspecs_unproc{31:38}],2);
% mean2014u = median([fspecs_unproc{39:46}],2);
% mean2015u = median([fspecs_unproc{47:73}],2);
% mean2016u = median([fspecs_unproc{74:79}],2);

% plot([mean2009u, mean2010u, mean2011u, mean2012u, mean2013u, mean2014u, mean2015u, mean2016u])
% legend('2009','2010','2011', '2012', '2013','2014','2015', '2016')

flat_intensity = zeros(1,length(fspecs_proc));
for i = 1:length(fspecs_proc)
    intensity = sum([fspecs_proc{i}]);
    flat_intensity(i) = intensity;
end


fdate = transpose([fdate{:}]);
%bad_flats = [3,6,21,24,28,29,30,51];
figure;
plot([fdate(:)], flat_intensity, '.-')
legend('flatfield intensity in sum of counts/s per pixel')

%fdate = transpose([fdate{:}]);

















% darkframes_copy = darkframes;
% dexptimes_copy = dexptimes;
% ddate_copy = ddate;
% dspecs_unproc_copy = dspecs_unproc;
% dspecs_proc_copy = dspecs_proc;
% 
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
% 
% 
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

