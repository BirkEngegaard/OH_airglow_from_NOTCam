function [useMaster, dlist, dexposureTime, nightMaster] = findDarks(file)
    [~, S] = getFitsData(file);
    nightMaster = 0;
    useMaster = 0;
    k = S.PrimaryData.Keywords;
    expMode = S.PrimaryData.Keywords{find(strcmp(k,'EXPMODE')), 2}; % Exposure mode e.g. 'frames 3.6 3'
    dexpMode = ['d', expMode];
    expMode_split = split(expMode, ' '); % split string into array for easier handling
    dexposureTime = str2double(expMode_split{2})*str2double(expMode_split{3});
    masterFile = importtxt('C:\Users\birke\Documents\prosjekt\ReduceSingleNOTfileBirk\NOTCamOverviewMaster.txt',',',2);
    [r,~] = size(masterFile);
    dlist = [];
  
    for i = 1:r
            name = masterFile{i,1};
            emode = masterFile{i,4};
            if isempty(name)
                continue;
            end
            if strcmp(name(1:6), file(1:6)) && contains(emode, 'dframes')
                dfile = masterFile{i,1}(1:10);
                [~, D] = getFitsData(dfile);
                kD = D.PrimaryData.Keywords;
                master_dexpMode = D.PrimaryData.Keywords{find(strcmp(k,'EXPMODE')), 2};
                if strcmp(master_dexpMode, dexpMode)
                    dlist = [dlist; dfile];
                end
            end
    end
  
    
    if isempty(dlist)
        nightMaster = 1;
        dexposureTime = [];
        k = 0;
        for i = 1:r
            name = masterFile{i,1};
            emode = masterFile{i,4};
            if isempty(name)
                continue;
            end
            if strcmp(name(1:6), file(1:6)) && contains(emode, 'dframes')
                k = k + 1;
                dfile = masterFile{i,1}(1:10);
                [~, D] = getFitsData(dfile);
                kD = D.PrimaryData.Keywords;
                master_dexpMode = D.PrimaryData.Keywords{find(strcmp(kD,'EXPMODE')), 2};
                master_dexpMode_split = split(master_dexpMode, ' '); % split string into array for easier handling
                dexposureTime_i = str2double(master_dexpMode_split{2})*str2double(master_dexpMode_split{3});
                dlist = [dlist; dfile];
                dexposureTime(k) = dexposureTime_i;
            end
            end
        
        
    end
    
    if isempty(dlist)
        useMaster = 1;
        nightMaster = 0;
    end
      
end
    
                