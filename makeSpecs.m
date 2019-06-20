%% Make specs for a given day
day = 'NCyi28';
files = findHfiles(day);
[rows, columns] = size(files);
for i = 1:rows
    SinglePrepare_dev(files(i, 1:10));
end
