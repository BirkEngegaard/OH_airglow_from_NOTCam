close all; clear all;
filter = 'H';
day = 'NCyi28';
path = ['D:/NOTArchive/',day, '/'];
%day = path(15:20);
dpath = path;
hpath = path;
[nil1, nil2, Hfilter, nil3] = masterfileanalysis('H', 'none');
for i = 1:size(Hfilter)
    if strcmp(Hfilter{i}(1:6), day)
        file = Hfilter{i};
        break
    end
end
%file = 'NCwb191106';
%day = file(1:6);

%% Make spec
data = fitsread(strcat(path,file,'.fits'),'Image');
S=fitsinfo(strcat(path,file,'.fits'));
k=S.PrimaryData.Keywords;

%% Get exposure time and choose the observation data point to be in the middle of observation
expMode = S.PrimaryData.Keywords{find(strcmp(k,'EXPMODE')), 2}; % Exposure mode e.g. 'frames 3.6 3'
dexpMode = ['d', expMode];
expMode = split(expMode, ' '); % split string into array for easier handling

try
    exposureTime = str2double(expMode{2})*str2double(expMode{3}); % Total time = time of ramp * num of ramps
catch
    warning('Couldn''t determine exposure time. Default value of 10.8 used')
    exposureTime = 10.8;
end

[flatfields, darkframes, Hfilter, master] = masterfileanalysis(filter, dexpMode);

hlist = [];
dlist = [];
for i = 1:size(flatfields)
    if strcmp(flatfields{i}(1:6), day)
        hlist = [hlist; flatfields{i}];
    end
end

for i = 1:size(darkframes)
    if strcmp(darkframes{i}(1:6), day)
        dlist = [dlist; darkframes{i}];
    end
end

[m, n] = size(dlist);
if m > 1
dlist(1,:) = [];
end

[dark,hot,flat,bad] = prepareBack(dpath,dlist,hpath,hlist);

data(data < -0) = -1;
data(abs(data) > (10000*exposureTime/10.8)) = -1; % Handle muon 10000

star = sum(data); %Handle Star
pos = find(star>500*(exposureTime/10.8)*1024);
if size(pos)>0
    starpos = pos(round(length(pos)/2));
    data(:,starpos-50:starpos+50) = -1;
end

data = data - dark;
data(bad == 1) = -1;
data(hot == 0) = -1;
data(abs(data) > 200*(exposureTime/10.8)) = -1; % Handle muon 300
bsumdata = data;
bsumdata(data>=0)=0;
bsumdata(data<0)=1;
baderased = sum(bsumdata,2);

data(data < 0) = 0;

%% Linearize

ldata = zeros(size(data)); %linearise
lflat = zeros(size(flat));
for i=1:1000
    for j=1:1000
        ldata(j,i)=data(round(f(i,j)),i); %Arcane magic, do not touch! it works! just believe me!! Kind regards Christoph from 2016
        lflat(j,i)=flat(round(f(i,j)),i);
    end  
end
telalt = k{find(strcmp(k,'TELALT')),2};
ldata=abs(ldata/sin(telalt));
spec = sum(ldata,2); %add up
spec = spec.*1024./(1024-baderased);
%spec = flipud(spec);


%% calibration
%M = dlmread(strcat(path,sprintf('%d-',53),file,'.txt'));
%spec=M(:,2);

pixel = linspace(1,1024,1024)';
pixRough = [896;862;822;810;785;771;590;551;540;502;367;323;307;288;277;251;239];
pix = zeros(size(pixRough));
i=1;
me = median(spec);
for p = pixRough'
    short = spec(p-10:p+10);
    y = linspace(p-10,p+10,21)';
    gauss = fittype('a*exp(-((x-m)/s)^2)+c','coeff',{'a','m','s','c'});
    f=fit(y,short,gauss,'StartPoint',[21000, p, 1.5, me],'Lower',[100,p-8,1,me*0.1],'Upper',[40000,p+8,4,4*me]);
    figure;
    plot(f,y,short);
    coeffs = coeffvalues(f);
    pix(i)=coeffs(2);
    i=i+1;
end
lam = [1.50555;1.51870;1.52876;1.53319;1.53951;1.54316;1.59725;1.60796;1.61281;1.62347;1.66923;1.68403;1.69033;1.69548;1.70082;1.70781;1.71229];

poly3 = fittype('a*x^3+b*x^2+c*x+d','coeff',{'a','b','c','d'});
f=fit(pix,lam,poly3,'StartPoint',[-1.6e-11,8.5e-08,3.4e-4,1.9]);
figure;
plot(f,pix,lam);

coeffs = coeffvalues(f);

p0 = coeffs(4);%1.944;
p1 = coeffs(3);%3.035e-4;
p2 = coeffs(2);%1.536e-7;
p3 = coeffs(1);%-5.705e-11;

%lambda = zeros(size(spec));
lambda = p0 + p1*pixel + p2*pixel.^2 + p3*pixel.^3;
% transpose(lambda);

residual = lam - (p0 + p1*pix + p2*pix.^2 + p3*pix.^3);
figure;
plot(pix,residual);
title('Residual');
xlabel('pixel number');
ylabel('\Delta\Lambda (\mum)');

figure;
plot(lambda,spec);
lambda = flipud(lambda);
save(strcat(path,'lambdaH.txt'),'lambda','-ascii');
disp([day, ' ', file, ' ', dexpMode])
%dlmwrite(strcat(path, sprintf('lambdaH2.txt')), lambda, '\n') % Birk's pref way of saving

%lambda1 = dlmread('C:/Users/chrifra/Desktop/BranchingRatio/prepareFits/lookup/lambdaKprime.txt');