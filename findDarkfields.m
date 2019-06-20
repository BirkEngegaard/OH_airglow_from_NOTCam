
function dlist = findDarkfields(path)
    path_to_file = [path, "Overview.txt"];
    overview = importtxt(path_to_file, ',', 2);
    dlist = [];
    dark1 = false;
    for line = 1:length(overview)
        overviewFile = convertStringsToChars(overview(line, 1));
        overviewExpMode = convertStringsToChars(overview(line, 3)); 
        
        if (contains(overviewExpMode, 'dfra'))
            if contains(overviewExpMode, ' 1') && (~dark1)
                dark1 = true;
                dlist = [dlist; overviewFile(1:10)];
            elseif (~(strcmp(overviewExpMode, 'Halogen 1')))
                dlist = [dlist; overviewFile(1:10)];
            end
        end
    end
    [r, c] = size(dlist);
    if (r > 1)
        dlist(1,:) = [];
    end
end
    