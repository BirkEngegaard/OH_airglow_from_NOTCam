function [q31counts, relerr31, q42counts, relerr42,q53counts,relerr53,q64counts,relerr64, SAT, specificDark, exposureTime, obsTime, airmass]...
    = getHfilterIntensities(file, masterDark, masterFlat,showFig)

    %% This function simply makes a call to SinglePrepare_dev.
    % Was too lazy to rename all the function calls for the H filter

    [q31counts, relerr31, q42counts, relerr42,q53counts,relerr53,q64counts,relerr64, SAT, specificDark, exposureTime, obsTime, airmass]...
    = SinglePrepare_dev(file, masterDark, masterFlat,showFig);

end