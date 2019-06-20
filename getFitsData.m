function [data, S] = getFitsData(file)
    %wd = 'C:\Users\birke\Documents\master\scripting\';
    wd = 'D:\NOTArchive\';
    path = [wd, file(1:6), '\'];
    data = fitsread(strcat(path,file,'.fits'),'Image');
    S=fitsinfo(strcat(path,file,'.fits'));
end