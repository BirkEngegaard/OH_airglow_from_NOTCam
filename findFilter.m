function filter = findFilter(overview, file)
    for line=1:size(overview)
        if (erase(overview{line}, '.fits') == file)
            filter = overview{line, 2};
            break
        end
    end
end