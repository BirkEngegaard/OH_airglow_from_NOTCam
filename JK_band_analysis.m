Kfilterdata = load('Kfilterdata.mat');
Ksfilterdata = load('Ksfilterdata.mat');
Jfilterdata = load('Jfilterdata.mat');

K_files = Kfilterdata.Kfilterdata;
Ks_files = Ksfilterdata.Kfilterdata;
J_files = Jfilterdata.Jfilterdata;

[r, ~] = size(K_files);
N_K = r;

[r, ~] = size(Ks_files);
N_Ks = r;

[r, ~] = size(J_files);
N_J = r;

useful_K_files = cell(1,N_K);
useful_K_files(1,:) = {0};
K_time = cell(1,N_K);
K_time(1,:) = {0};

useful_Ks_files = cell(1,N_Ks);
useful_Ks_files(1,:) = {0};
Ks_time = cell(1,N_Ks);
Ks_time(1,:) = {0};

useful_J_files = cell(1,N_J);
useful_J_files(1,:) = {0};
J_time = cell(1,N_J);
J_time(1,:) = {0};

%% K filter
k_K = 0;
f = waitbar(0,'Scanning K filter files');
for i = 1:N_K
    %fault = 0;
    [~, S] = getFitsData(K_files{i});
    k=S.PrimaryData.Keywords;
    [obsindex,~] = find(strcmp(k,'OBS_MODE'));
    [aindex1, ~] = find(strcmp(k,'APERTUR'));
    [aindex2, ~] = find(strcmp(k,'NCAPRNM'));
    [findex,~] = find(strcmp(k,'NCFLTNM2'));
    [findex2,~] = find(strcmp(k,'FILT2'));
    if ~strcmp(k{obsindex,2}, 'Spectroscopy')
        continue;
    elseif ~isempty(aindex1)
        if ~strcmp(k{aindex1,2}, '128mu WF slit')
            continue;
        end
    elseif ~isempty(aindex2)
        if ~strcmp(k{aindex2,2}, '128mu WF slit')
            continue;
        end
    end
    k_K = k_K + 1;
    useful_K_files{k_K} = K_files{i};
    [oindex, ~] = find(strcmp(S.PrimaryData.Keywords,'DATE-OBS'));
    obsTime = datetime(replace(S.PrimaryData.Keywords{oindex, 2}, 'T', ' '), 'Format', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
    K_time{k_K} = obsTime;
    waitbar(i/N_K,f, '(1/3) Scanning K filter files')
end
close(f)

%% Ks filter
k_Ks = 0;
f = waitbar(0,'(2/3) Scanning Ks filter files');
for i = 1:N_Ks
    %fault = 0;
    [~, S] = getFitsData(Ks_files{i});
    k=S.PrimaryData.Keywords;
    [obsindex,~] = find(strcmp(k,'OBS_MODE'));
    [aindex1, ~] = find(strcmp(k,'APERTUR'));
    [aindex2, ~] = find(strcmp(k,'NCAPRNM'));
    [findex,~] = find(strcmp(k,'NCFLTNM2'));
    [findex2,~] = find(strcmp(k,'FILT2'));
    if ~strcmp(k{obsindex,2}, 'Spectroscopy')
        continue;
    elseif ~isempty(aindex1)
        if ~strcmp(k{aindex1,2}, '128mu WF slit')
            continue;
        end
    elseif ~isempty(aindex2)
        if ~strcmp(k{aindex2,2}, '128mu WF slit')
            continue;
        end
    end
    k_Ks = k_Ks + 1;
    useful_Ks_files{k_Ks} = Ks_files{i};
    [oindex, ~] = find(strcmp(S.PrimaryData.Keywords,'DATE-OBS'));
    obsTime = datetime(replace(S.PrimaryData.Keywords{oindex, 2}, 'T', ' '), 'Format', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
    Ks_time{k_Ks} = obsTime;
    waitbar(i/N_Ks,f, '(2/3) Scanning Ks filter files')
end


close(f)

%% J filter
k_J = 0;
f = waitbar(0,'(3/3) Scanning J filter files');
for i = 1:N_J
    %fault = 0;
    [~, S] = getFitsData(J_files{i});
    k=S.PrimaryData.Keywords;
    [obsindex,~] = find(strcmp(k,'OBS_MODE'));
    [aindex1, ~] = find(strcmp(k,'APERTUR'));
    [aindex2, ~] = find(strcmp(k,'NCAPRNM'));
    [findex,~] = find(strcmp(k,'NCFLTNM2'));
    [findex2,~] = find(strcmp(k,'FILT2'));
    if ~strcmp(k{obsindex,2}, 'Spectroscopy')
        continue;
    elseif ~isempty(aindex1)
        if ~strcmp(k{aindex1,2}, '128mu WF slit')
            continue;
        end
    elseif ~isempty(aindex2)
        if ~strcmp(k{aindex2,2}, '128mu WF slit')
            continue;
        end
    end
    k_J = k_J + 1;
    useful_Ks_files{k_J} = J_files{i};
    [oindex, ~] = find(strcmp(S.PrimaryData.Keywords,'DATE-OBS'));
    obsTime = datetime(replace(S.PrimaryData.Keywords{oindex, 2}, 'T', ' '), 'Format', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
    J_time{k_J} = obsTime;
    waitbar(i/N_J,f, '(3/3) Scanning J filter files')
end


close(f)

useful_K_files = useful_K_files(1:k_K);
K_time = K_time(1:k_K);

useful_Ks_files = useful_Ks_files(1:k_Ks);
Ks_time = Ks_time(1:k_Ks);

useful_J_files = useful_J_files(1:k_J);
J_time = J_time(1:k_J);

figure;
    histogram2([day([K_time{:}], 'dayofyear')], [year([K_time{:}])],[365,9], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year', 'FontSize', 14)
    ylabel('Year', 'FontSize', 14)
    title('K filter')
    colorbar
figure;
    histogram2([day([Ks_time{:}], 'dayofyear')], [year([Ks_time{:}])],[365,9], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year')
    ylabel('Year')
    title('Ks filter')
    colorbar
    
 figure;
    histogram2([day([J_time{:}], 'dayofyear')], [year([J_time{:}])],[365,8], 'ShowEmptyBins','off', 'DisplayStyle', 'tile')
    xlabel('Day of year', 'FontSize', 14)
    ylabel('Year', 'FontSize', 14)
    title('J filter')
    colorbar


