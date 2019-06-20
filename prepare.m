function [q31counts,relerr31,q42counts,relerr42,q53counts,relerr53,q64counts,relerr64,...
    SAT, specificDark, exposureTime, obsTime, airmass, rampMode]...
    = prepare(path,file,dpath,dlist,hpath,hlist,trans,filter, figpath, S, data, masterDark, masterFlat,showFig)

%[dark,hot,flat,bad] = prepareBack(dpath,dlist,hpath,hlist,file);
[dark,hot,flat,bad, useMasterDark, nightMaster] = prepareBack_master(dpath,file,filter, masterDark, masterFlat);

[q31counts, relerr31, q42counts, relerr42,q53counts,relerr53,q64counts,relerr64,...
    SAT, specificDark, exposureTime, obsTime, airmass, rampMode]...
    = preparePic(path,file,trans,filter,dark,hot,flat,bad,figpath,S,data,showFig, useMasterDark, nightMaster);

% try
% preparePic(path,file,trans,filter,dark,hot,flat,bad,figpath, S, data);
% catch
% errstr = sprintf('Error in preparePic for file: %s', file);
% error(errstr)
% end
end