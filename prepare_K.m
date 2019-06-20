function [q97counts, relerr97, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = prepare_K(path,file,filter, figpath, S, data, masterDark, masterFlat,showFig)

%[dark,hot,flat,bad] = prepareBack(dpath,dlist,hpath,hlist,file);
[dark,hot,flat,bad, useMasterDark, nightMaster] = prepareBack_master(path,file,filter, masterDark, masterFlat);

[q97counts, relerr97, specificDark, exposureTime, obsTime, airmass, SAT, rampMode]...
    = preparePic_K(path,file,filter,dark,hot,flat,bad,figpath,S,data,showFig, useMasterDark, nightMaster);
end