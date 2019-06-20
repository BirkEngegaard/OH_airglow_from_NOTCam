function [flatfields, darkframes, Jfilterdata, master] = masterfileanalysisJ(filter, dexpmode)
filter = 'J';
dexpmode = 'none';
master = importtxt('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NOTCamOverviewMaster.txt',',',2);

%fileID = fopen('H_filter_data.txt','w');
Jfilterdata = {};
k = 1;
darkframes = {};
d = 1;

flatfields = {};
f = 1;

if strcmp(dexpmode,'none')
   for line = 1:size(master)
    if strcmp(master{line,2},filter) && contains(master{line,4}, 'frames ') && ~ contains(master{line,4}, 'dframes ')
        file = master{line,1}(1:10);
        Jfilterdata{k} = file;
        k = k + 1;
    end
   end
Jfilterdata = transpose(Jfilterdata);
else

for line = 1:size(master)
    if strcmp(master{line,2},filter) && contains(master{line,4}, 'frames ') && ~ contains(master{line,4}, 'dframes ')
        file = master{line,1}(1:10);
        Jfilterdata{h} = file;
        h = h + 1;
    end
    if contains(master{line,4}, dexpmode)
        dfile = master{line,1}(1:10);
        darkframes{d} = dfile;
        d = d + 1;
    end
    if contains(master{line,3}, 'Halogen') && ~strcmp(master{line,3}, 'Halogen 1') && strcmp(master{line,2},filter)
        hfile = master{line,1}(1:10);
        flatfields{f} = hfile;
        f = f + 1;
    end
end

Jfilterdata = transpose(Jfilterdata);
flatfields = transpose(flatfields);
darkframes = transpose(darkframes);

flat_days = {};
j = 1;
for i = 1:size(flatfields)
    if i == 1
        flat_days{j} = flatfields{i}(1:6);
    elseif ~strcmp(flatfields{i}(1:6), flatfields{i-1}(1:6))
        j = j + 1;
        flat_days{j} = flatfields{i}(1:6);
    end
end
end
end
%Hfilterdata = textscan(fID,'%s');
        