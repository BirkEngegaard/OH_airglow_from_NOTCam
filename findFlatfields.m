
function hlist = findFlatfields(path, filter)
    path_to_file = [path, "Overview.txt"];
    overview = importtxt(path_to_file, ',', 2);
    hlist = [];
    halogen1 = false;
    for line = 1:length(overview)
        overviewFile = convertStringsToChars(overview(line, 1));
        overviewFilter = convertStringsToChars(overview(line, 2));
        overviewExpMode = convertStringsToChars(overview(line, 3)); 
        if (strcmp(overviewFilter, filter)) && (contains(overviewExpMode, 'Halogen'))
            if strcmp(overviewExpMode, 'Halogen 1') && (~halogen1)
                halogen1 = true;
                hlist = [hlist; overviewFile(1:10)];
            elseif (~(strcmp(overviewExpMode, 'Halogen 1')))
                hlist = [hlist; overviewFile(1:10)];
            end
        end
    end
    [r, c] = size(hlist);
    if (r > 1)
        hlist(1,:) = [];
    end
end
    