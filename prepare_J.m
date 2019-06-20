function [q74counts, relerr74, q85counts, relerr85, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = prepare_J(path,file,filter, figpath, S, data, masterDark, masterFlat,showFig)

%[dark,hot,flat,bad] = prepareBack(dpath,dlist,hpath,hlist,file);
[dark,hot,flat,bad, useMasterDark, nightMaster] = prepareBack_master(path,file,filter, masterDark, masterFlat);

[q74counts, relerr74, q85counts, relerr85, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = preparePic_J(path,file,filter,dark,hot,flat,bad,figpath,S,data,showFig, useMasterDark, nightMaster);
end