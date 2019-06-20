H_month_dist = zeros(1,12);
K_month_dist = zeros(1,12);

for i = 1:12
    H_month_dist(i) = sum(month([H_time{:}]) == i);
    K_month_dist(i) = sum(month([K_time{:}]) == i);
end

figure;
histogram(month([H_time{:}]), 0.5:11.5, 'FaceColor','r')
hold on;
histogram(month([K_time{:}]), 0.5:11.5,'FaceColor', 'b')
xlim([0.5,12.5])
xticks(1:12)
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
legend('$H$ band', '$K$ band')
xlabel('Month')
ylabel('Number of observations')

figure
bar([0.85:11.85], H_month_dist, 0.30)
hold on;
bar([1.15:12.15], K_month_dist, 0.30)
xlim([0.70,12.3])
xticks(1:12)
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
legend('$H$ band', '$K$ band')
xlabel('Month')
ylabel('Number of observations')
ylim([0,150])

