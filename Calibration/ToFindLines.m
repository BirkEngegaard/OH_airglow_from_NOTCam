clear all; close all;

trans = 86;

PHYS_kB = 0.6950356; % in cm^-1/K
Mpar = dlmread(sprintf('C:/Users/chrifra/Desktop/Papers/Report/Prog/Data/FitData/%d.par',trans));
%l=1./Mpar(:,1)*10000;
l=Mpar(:,1);
Eu=Mpar(:,2);
A=Mpar(:,3);
gu=Mpar(:,4);
I=gu.*exp(-1.*Eu/(PHYS_kB.*200)).*A.*1e25*1e48;
pix = 1:size(l);

figure;
plot(pix,I);
axis([150,250,0,max(I)]);