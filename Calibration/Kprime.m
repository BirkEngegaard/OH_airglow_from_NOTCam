close all; clear all;

path = 'D:/NOTArchive/NCzf20/';
file = 'NCzf200400';

M = dlmread(strcat(path,sprintf('%d-',86),file,'.txt'));
spec=M(:,2);
yes = max(spec)*2;
pos = zeros(size(spec));
for i=20:size(spec)-20
    if spec(i)==max(spec(i-10:i+10))
        pos(i)=1;
    end
end

linesfound = sum(pos); 

figure;
pixel = linspace(1,1024,1024)';
plot(pixel,spec); hold on;
plot(pixel,pos*yes);
axis([1,1024,0,max(spec(200:800))]);

pixRough = [112;177;249;284;324;366;557;626;662;701;744];
pix = zeros(size(pixRough));
i=1;
me = median(spec);
for p = pixRough'
    short = spec(p-10:p+10);
    y = linspace(p-10,p+10,21)';
    gauss = fittype('a*exp(-((x-m)/s)^2)+c','coeff',{'a','m','s','c'});
    f=fit(y,short,gauss,'StartPoint',[21000, p, 1.5, me],'Lower',[500,p-10,1,me*0.1],'Upper',[40000,p+10,4,4*me])
    figure;
    plot(f,y,short);
    coeffs = coeffvalues(f);
    pix(i)=coeffs(2);
    i=i+1;
end
lam = [1.97720;2.000814;2.02759;2.04127;2.05635;2.07291;2.15073;2.18023;2.19556;2.21255;2.23127];

poly3 = fittype('a*x^3+b*x^2+c*x+d','coeff',{'a','b','c','d'});
f=fit(pix,lam,poly3,'StartPoint',[-1.6e-11,8.5e-08,3.4e-4,1.9]);
plot(f,pix,lam);

coeffs = coeffvalues(f);

p0 = coeffs(4);%1.944;
p1 = coeffs(3);%3.035e-4;
p2 = coeffs(2);%1.536e-7;
p3 = coeffs(1);%-5.705e-11;

%lambda = zeros(size(spec));
lambda = p0 + p1*pixel + p2*pixel.^2 + p3*pixel.^3;
%transpose(lambda);

residual = lam - (p0 + p1*pix + p2*pix.^2 + p3*pix.^3);
figure;
plot(pix,residual);
title('Residual');
xlabel('pixel number');
ylabel('\Delta\Lambda (\mum)');

figure;
plot(lambda,spec);
% save(strcat(path,'lambdaKprime.txt'),'lambda','-ascii');

%lambda1 = dlmread('C:/Users/chrifra/Desktop/BranchingRatio/prepareFits/lookup/lambdaKprime.txt');