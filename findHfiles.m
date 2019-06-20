function Hfile = findHfiles(day)
    master = importtxt('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NOTCamOverviewMaster.txt',',',2);
    filter = 'H';
    Hfilterdata = {};
    h = 1;  
    for line = 1:size(master)
        if strcmp(master{line,2},'H') && contains(master{line,4}, 'frames ') && ~ contains(master{line,4}, 'dframes ')
            file = master{line,1}(1:10);
            Hfilterdata{h} = file;
            h = h + 1;
        end
    end
    
    j = 0;
    Hfile = {};
    Hfilterdata = transpose(Hfilterdata);
    for i = 1:size(Hfilterdata)
        if strcmp(Hfilterdata{i}(1:6), day)
            j = j + 1;
            Hfile{j,1} = Hfilterdata{i};
        end
    end
end