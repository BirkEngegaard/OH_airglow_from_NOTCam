load('OHintensity_data_HK_280319.mat')

intensity31_year_iteration = intensity31_year;
intensity42_year_iteration = intensity42_year;
intensity53_year_iteration = intensity53_year;
intensity64_year_iteration = intensity64_year;
intensity97_year_iteration = intensity97_year;

% intensity31_seasonRemoved = intensity31;
% intensity42_seasonRemoved = intensity42;
% intensity53_seasonRemoved = intensity53;
% intensity64_seasonRemoved = intensity64;
% intensity97_seasonRemoved = intensity97;

for iteration = 1:1

%% In terms of month -> change of variable to noy

intensity31_moy_expanded = [intensity31_year_iteration intensity31_year_iteration intensity31_year_iteration];
intensity42_moy_expanded = [intensity42_year_iteration intensity42_year_iteration intensity42_year_iteration];
intensity53_moy_expanded = [intensity53_year_iteration intensity53_year_iteration intensity53_year_iteration];
intensity64_moy_expanded = [intensity64_year_iteration intensity64_year_iteration intensity64_year_iteration];
intensity97_moy_expanded = [intensity97_year_iteration intensity97_year_iteration intensity97_year_iteration];
moy_expanded = [0.5:35.5];

t = ~isnan(intensity31_moy_expanded);
t_K = ~isnan(intensity97_moy_expanded);
noy_moy = linspace(0,12,365);
noy_expanded = [noy_moy noy_moy+12 noy_moy+12+12];

intensity31_moy_spline_expanded = spline(moy_expanded(t), intensity31_moy_expanded(t), noy_expanded);
intensity42_moy_spline_expanded = spline(moy_expanded(t), intensity42_moy_expanded(t), noy_expanded);
intensity53_moy_spline_expanded = spline(moy_expanded(t), intensity53_moy_expanded(t), noy_expanded);
intensity64_moy_spline_expanded = spline(moy_expanded(t), intensity64_moy_expanded(t), noy_expanded);
intensity97_moy_spline_expanded = spline(moy_expanded(t_K), intensity97_moy_expanded(t_K), noy_expanded);
intensity31_noy_spline = intensity31_moy_spline_expanded(366:730); % Intepret each monthly value as the value halfway into the month
intensity42_noy_spline = intensity42_moy_spline_expanded(366:730);
intensity53_noy_spline = intensity53_moy_spline_expanded(366:730);
intensity64_noy_spline = intensity64_moy_spline_expanded(366:730);
intensity97_noy_spline = intensity97_moy_spline_expanded(366:730);

figure;
plot(0.5:11.5, intensity31_year_iteration, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
hold on
plot(0.5:11.5, intensity42_year_iteration, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(0.5:11.5, intensity53_year_iteration, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(0.5:11.5, intensity64_year_iteration, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(0.5:11.5, intensity97_year_iteration, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(noy_moy, intensity31_noy_spline, 'LineWidth', 3)
plot(noy_moy, intensity42_noy_spline, 'LineWidth', 3)
plot(noy_moy, intensity53_noy_spline, 'LineWidth', 3)
plot(noy_moy, intensity64_noy_spline, 'LineWidth', 3)
plot(noy_moy, intensity97_noy_spline, 'LineWidth', 3)
xticks(0.5:11.5)
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
legend('Q$(3,1)$','Q$(4,2)$','Q$(5,3)$','Q$(6,4)$','Q$(9,7)$','Q$(3,1)$ spline','Q$(4,2)$ spline','Q$(5,3)$ spline','Q$(6,4)$ spline','Q$(9,7)$ spline')
xlabel('Month')
ylabel('Intensity [a.u.]')
title('Spline interp. of seasonal cycle')

figure;
plot(moy_expanded,intensity31_moy_expanded, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
hold on
plot(moy_expanded,intensity42_moy_expanded, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(moy_expanded,intensity53_moy_expanded, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(moy_expanded,intensity64_moy_expanded, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(moy_expanded,intensity97_moy_expanded, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(noy_expanded, intensity31_moy_spline_expanded, 'LineWidth', 3)
plot(noy_expanded,intensity42_moy_spline_expanded, 'LineWidth', 3)
plot(noy_expanded,intensity53_moy_spline_expanded, 'LineWidth', 3)
plot(noy_expanded,intensity64_moy_spline_expanded, 'LineWidth', 3)
plot(noy_expanded, intensity97_moy_spline_expanded, 'LineWidth', 3)
%xticks(0.5:11.5)
%xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
legend('raw','raw','raw','raw','raw','31','42','53','64','97')
%xlabel('Month')
ylabel('Intensity [a.u.]')
title('Spline interp. of seasonal cycle')

%% Make new data point with season removed

N = length(intensity31);
N_K = length(intensity97);
intensity31_seasonRemoved = zeros(1,N);
intensity42_seasonRemoved = zeros(1,N);
intensity53_seasonRemoved = zeros(1,N);
intensity64_seasonRemoved = zeros(1,N);
intensity97_seasonRemoved = zeros(1,N_K);

for i = 1:N
    temp_noy = H_noy(i);
    intensity31_seasonRemoved(i) = intensity31(i) - intensity31_noy_spline(temp_noy);
    intensity42_seasonRemoved(i) = intensity42(i) - intensity42_noy_spline(temp_noy);
    intensity53_seasonRemoved(i) = intensity53(i) - intensity53_noy_spline(temp_noy);
    intensity64_seasonRemoved(i) = intensity64(i) - intensity64_noy_spline(temp_noy);
end

for i = 1:N_K
    temp_noy = K_noy(i);
    intensity97_seasonRemoved(i) = intensity97(i) - intensity97_noy_spline(temp_noy);
end

% for i = 1:N
%     temp_noy = H_noy(i);
%     intensity31_seasonRemoved(i) = intensity31_seasonRemoved(i) - intensity31_noy_spline(temp_noy);
%     intensity42_seasonRemoved(i) = intensity42_seasonRemoved(i) - intensity42_noy_spline(temp_noy);
%     intensity53_seasonRemoved(i) = intensity53_seasonRemoved(i) - intensity53_noy_spline(temp_noy);
%     intensity64_seasonRemoved(i) = intensity64_seasonRemoved(i) - intensity64_noy_spline(temp_noy);
% end
% 
% for i = 1:N_K
%     temp_noy = K_noy(i);
%     intensity97_seasonRemoved(i) = intensity97_seasonRemoved(i) - intensity97_noy_spline(temp_noy);
% end

lst_hours_lowres = [17 18 19 20 21 22 23 0 1 2 3 4 5 6];
intensity31_lst_lowres_seasonRemoved = zeros(1,14);
intensity42_lst_lowres_seasonRemoved = zeros(1,14);
intensity53_lst_lowres_seasonRemoved = zeros(1,14);
intensity64_lst_lowres_seasonRemoved = zeros(1,14);
intensity97_lst_lowres_seasonRemoved = zeros(1,14);
err31_lst_lowres_seasonRemoved = zeros(1,14);
err42_lst_lowres_seasonRemoved = zeros(1,14);
err53_lst_lowres_seasonRemoved = zeros(1,14);
err64_lst_lowres_seasonRemoved = zeros(1,14);
err97_lst_lowres_seasonRemoved = zeros(1,14);

for i = 1:14
    temp_int31 = zeros(1, length(H_LST));
    temp_int42 = zeros(1, length(H_LST));
    temp_int53 = zeros(1, length(H_LST));
    temp_int64 = zeros(1, length(H_LST));
    k = 0;
    for j = 1:length(H_LST)
        if str2double(H_LST{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int31(k) = intensity31_seasonRemoved(j);
            temp_int42(k) = intensity42_seasonRemoved(j);
            temp_int53(k) = intensity53_seasonRemoved(j);
            temp_int64(k) = intensity64_seasonRemoved(j);
        end
    end
    temp_int31 = temp_int31(1:k);
    temp_int42 = temp_int42(1:k);
    temp_int53 = temp_int53(1:k);
    temp_int64 = temp_int64(1:k);
    intensity31_lst_lowres_seasonRemoved(i) = mean(temp_int31);
    intensity42_lst_lowres_seasonRemoved(i) = mean(temp_int42);
    intensity53_lst_lowres_seasonRemoved(i) = mean(temp_int53);
    intensity64_lst_lowres_seasonRemoved(i) = mean(temp_int64);
    err31_lst_lowres_seasonRemoved(i) = std(temp_int31)/sqrt(length(temp_int31));
    err42_lst_lowres_seasonRemoved(i) = std(temp_int42)/sqrt(length(temp_int42));
    err53_lst_lowres_seasonRemoved(i) = std(temp_int53)/sqrt(length(temp_int53));
    err64_lst_lowres_seasonRemoved(i) = std(temp_int64)/sqrt(length(temp_int64));
    
    temp_int97 = zeros(1, length(K_LST));
    k = 0;
    for j = 1:length(K_LST)
        if str2double(K_LST{j}(12:13)) == lst_hours_lowres(i)
            k = k + 1;
            temp_int97(k) = intensity97_seasonRemoved(j);
        end
    end
    temp_int97 = temp_int97(1:k);
    intensity97_lst_lowres_seasonRemoved(i) = mean(temp_int97);
    err97_lst_lowres_seasonRemoved(i) = std(temp_int97)/sqrt(length(temp_int97));
end
    


%% Remove night from season iteration
intensity31_lst_lowres_seasonRemoved_expanded =...
    [intensity31_lst_lowres_seasonRemoved...
    intensity31_lst_lowres_seasonRemoved intensity31_lst_lowres_seasonRemoved];
intensity42_lst_lowres_seasonRemoved_expanded =...
    [intensity42_lst_lowres_seasonRemoved...
    intensity42_lst_lowres_seasonRemoved intensity42_lst_lowres_seasonRemoved];
intensity53_lst_lowres_seasonRemoved_expanded =...
    [intensity53_lst_lowres_seasonRemoved...
    intensity53_lst_lowres_seasonRemoved intensity53_lst_lowres_seasonRemoved];
intensity64_lst_lowres_seasonRemoved_expanded =...
    [intensity64_lst_lowres_seasonRemoved...
    intensity64_lst_lowres_seasonRemoved intensity64_lst_lowres_seasonRemoved];
intensity97_lst_lowres_seasonRemoved_expanded =...
    [intensity97_lst_lowres_seasonRemoved...
    intensity97_lst_lowres_seasonRemoved intensity97_lst_lowres_seasonRemoved];

hp17lst = [0.5:13.5];
hp17lst_expanded = [hp17lst hp17lst+14 hp17lst+28];

t = ~isnan(intensity31_lst_lowres_seasonRemoved_expanded);
t_K = ~isnan(intensity97_lst_lowres_seasonRemoved_expanded);

xq = linspace(0, 42, 2520); % one index for each minute after 17 LST

intensity31_lst_lowres_seasonRemoved_expanded_spline = spline(hp17lst_expanded(t), intensity31_lst_lowres_seasonRemoved_expanded(t), xq);
intensity42_lst_lowres_seasonRemoved_expanded_spline = spline(hp17lst_expanded(t), intensity42_lst_lowres_seasonRemoved_expanded(t), xq);
intensity53_lst_lowres_seasonRemoved_expanded_spline = spline(hp17lst_expanded(t), intensity53_lst_lowres_seasonRemoved_expanded(t), xq);
intensity64_lst_lowres_seasonRemoved_expanded_spline = spline(hp17lst_expanded(t), intensity64_lst_lowres_seasonRemoved_expanded(t), xq);
intensity97_lst_lowres_seasonRemoved_expanded_spline = spline(hp17lst_expanded(t_K), intensity97_lst_lowres_seasonRemoved_expanded(t_K), xq);
intensity31_lst_lowres_seasonRemoved_spline = intensity31_lst_lowres_seasonRemoved_expanded_spline(841:1680);
intensity42_lst_lowres_seasonRemoved_spline = intensity42_lst_lowres_seasonRemoved_expanded_spline(841:1680);
intensity53_lst_lowres_seasonRemoved_spline = intensity53_lst_lowres_seasonRemoved_expanded_spline(841:1680);
intensity64_lst_lowres_seasonRemoved_spline = intensity64_lst_lowres_seasonRemoved_expanded_spline(841:1680);
intensity97_lst_lowres_seasonRemoved_spline = intensity97_lst_lowres_seasonRemoved_expanded_spline(841:1680);

figure;
plot(hp17lst, intensity31_lst_lowres_seasonRemoved, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
hold on;
plot(hp17lst, intensity42_lst_lowres_seasonRemoved, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(hp17lst, intensity53_lst_lowres_seasonRemoved, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(hp17lst, intensity64_lst_lowres_seasonRemoved, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(hp17lst, intensity97_lst_lowres_seasonRemoved, 'b:', 'Marker', 'Hexagram', 'MarkerSize', 10)
plot(linspace(0,14,840),intensity31_lst_lowres_seasonRemoved_spline, 'LineWidth', 3)
plot(linspace(0,14,840),intensity42_lst_lowres_seasonRemoved_spline, 'LineWidth', 3)
plot(linspace(0,14,840),intensity53_lst_lowres_seasonRemoved_spline, 'LineWidth', 3)
plot(linspace(0,14,840),intensity64_lst_lowres_seasonRemoved_spline, 'LineWidth', 3)
plot(linspace(0,14,840),intensity97_lst_lowres_seasonRemoved_spline, 'LineWidth', 3)
xticks(0.5:13.5)
xticklabels(lst_hours_lowres)
legend('Q$(3,1)$','Q$(4,2)$','Q$(5,3)$','Q$(6,4)$','Q$(9,7)$','Q$(3,1)$ spline','Q$(4,2)$ spline','Q$(5,3)$ spline','Q$(6,4)$ spline','Q$(9,7)$ spline')
title('Spline fit to LST data. Season removed')


hplst17s = linspace(0,14,840);

%% Make new data points with night removed
intensity31_nightRemoved = zeros(1, length(intensity31));
intensity42_nightRemoved = zeros(1, length(intensity42));
intensity53_nightRemoved = zeros(1, length(intensity53));
intensity64_nightRemoved = zeros(1, length(intensity64));
intensity97_nightRemoved = zeros(1, length(intensity97));

%H filter
for i = 1:length(H_LST)
    temp_hour = str2double(H_LST{i}(12:13));
    temp_minute = str2double(H_LST{i}(15:16));
    
    temp_hp17lst = temp_hour - 17;
    if temp_hp17lst < 0
        temp_hp17lst = 7 + temp_hour;
    end

    hp17lst_index = 60*temp_hp17lst + temp_minute;
    
    intensity31_nightRemoved(i) =...
        intensity31(i) - intensity31_lst_lowres_seasonRemoved_spline(hp17lst_index);
    intensity42_nightRemoved(i) =...
        intensity42(i) - intensity42_lst_lowres_seasonRemoved_spline(hp17lst_index);
    intensity53_nightRemoved(i) =...
        intensity53(i) - intensity53_lst_lowres_seasonRemoved_spline(hp17lst_index);
    intensity64_nightRemoved(i) =...
        intensity64(i) - intensity64_lst_lowres_seasonRemoved_spline(hp17lst_index);
end
%K filter
for i = 1:length(K_LST)
    temp_hour = str2double(K_LST{i}(12:13));
    temp_minute = str2double(K_LST{i}(15:16));
    
    temp_hp17lst = temp_hour - 17;
    if temp_hp17lst < 0
        temp_hp17lst = 7 + temp_hour;
    end

    hp17lst_index = 60*temp_hp17lst + temp_minute;
    
    intensity97_nightRemoved(i) =...
        intensity97(i) - intensity97_lst_lowres_seasonRemoved_spline(hp17lst_index);
end



%% Seasonal cycle, night removed
intensity31_year_nightRemoved = zeros(1,12);
intensity42_year_nightRemoved = zeros(1,12);
intensity53_year_nightRemoved = zeros(1,12);
intensity64_year_nightRemoved = zeros(1,12);
intensity97_year_nightRemoved = zeros(1,12);
err31_year_nightRemoved = zeros(1,12);
err42_year_nightRemoved = zeros(1,12);
err53_year_nightRemoved = zeros(1,12);
err64_year_nightRemoved = zeros(1,12);
err97_year_nightRemoved = zeros(1,12);

for i = 1:12
    intensity31_year_nightRemoved(i) = mean(intensity31_nightRemoved(month([time{:}]) == i));
    intensity42_year_nightRemoved(i) = mean(intensity42_nightRemoved(month([time{:}]) == i));
    intensity53_year_nightRemoved(i) = mean(intensity53_nightRemoved(month([time{:}]) == i));
    intensity64_year_nightRemoved(i) = mean(intensity64_nightRemoved(month([time{:}]) == i));
    intensity97_year_nightRemoved(i) = mean(intensity97_nightRemoved(month([K_time{:}]) == i));
    
    err31_year_nightRemoved(i) = std(intensity31_nightRemoved(month([time{:}]) == i))/...
        sqrt(length(intensity31_nightRemoved(month([time{:}]) == i)));
    err42_year_nightRemoved(i) = std(intensity42_nightRemoved(month([time{:}]) == i))/...
        sqrt(length(intensity42_nightRemoved(month([time{:}]) == i)));
    err53_year_nightRemoved(i) = std(intensity53_nightRemoved(month([time{:}]) == i))/...
        sqrt(length(intensity53_nightRemoved(month([time{:}]) == i)));
    err64_year_nightRemoved(i) = std(intensity64_nightRemoved(month([time{:}]) == i))/...
        sqrt(length(intensity64_nightRemoved(month([time{:}]) == i)));
    err97_year_nightRemoved(i) = std(intensity97_nightRemoved(month([time{:}]) == i))/...
        sqrt(length(intensity97_nightRemoved(month([time{:}]) == i)));
end
return
if iteration == 1
    figure
    hold on;
    errorbar(1:14, intensity31_lst_lowres_seasonRemoved, err31_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity42_lst_lowres_seasonRemoved, err42_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity53_lst_lowres_seasonRemoved, err53_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity64_lst_lowres_seasonRemoved, err64_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity97_lst_lowres_seasonRemoved, err97_lst_lowres_seasonRemoved)
    xticks(1:14);
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    title('Nightly variation when season is removed')
    legend('31','42','53','64','97')

    figure;
    hold on;
    errorbar(0.5:11.5, intensity31_year_nightRemoved, err31_year_nightRemoved)
    errorbar(0.5:11.5, intensity42_year_nightRemoved, err42_year_nightRemoved)
    errorbar(0.5:11.5, intensity53_year_nightRemoved, err53_year_nightRemoved)
    errorbar(0.5:11.5, intensity64_year_nightRemoved, err64_year_nightRemoved)
    errorbar(0.5:11.5, intensity97_year_nightRemoved, err97_year_nightRemoved)
    xticks(0.5:11.5)
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    ylabel('Intensity [a.u.]')
    title('Annual variation, nightgly variation removed')
    legend('31','42','53','64','97')
end



intensity31_year_iteration = intensity31_year_nightRemoved;
intensity42_year_iteration = intensity42_year_nightRemoved;
intensity53_year_iteration = intensity53_year_nightRemoved;
intensity64_year_iteration = intensity64_year_nightRemoved;
intensity97_year_iteration = intensity97_year_nightRemoved;

end

if iteration ~= 1
figure
    hold on;
    errorbar(1:14, intensity31_lst_lowres_seasonRemoved, err31_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity42_lst_lowres_seasonRemoved, err42_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity53_lst_lowres_seasonRemoved, err53_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity64_lst_lowres_seasonRemoved, err64_lst_lowres_seasonRemoved)
    errorbar(1:14, intensity97_lst_lowres_seasonRemoved, err97_lst_lowres_seasonRemoved)
    xticks(1:14);
    xticklabels(lst_hours_lowres)
    xlabel('Local solar time')
    ylabel('Intensity [a.u.]')
    title('Nightly variation when season is removed')
    legend('31','42','53','64','97')

    figure;
    hold on;
    errorbar(0.5:11.5, intensity31_year_nightRemoved, err31_year_nightRemoved)
    errorbar(0.5:11.5, intensity42_year_nightRemoved, err42_year_nightRemoved)
    errorbar(0.5:11.5, intensity53_year_nightRemoved, err53_year_nightRemoved)
    errorbar(0.5:11.5, intensity64_year_nightRemoved, err64_year_nightRemoved)
    errorbar(0.5:11.5, intensity97_year_nightRemoved, err97_year_nightRemoved)
    xticks(0.5:11.5)
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    ylabel('Intensity [a.u.]')
    title('Annual variation, nightgly variation removed')
    legend('31','42','53','64','97')
end