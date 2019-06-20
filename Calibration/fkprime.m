function lambda = fkprime(spec)

pixRough = [284;324;366;557;626;662;701;744];
pix = zeros(size(pixRough));
i=1;
me = median(spec);
for p = pixRough'
    short = spec(p-10:p+10);
    y = linspace(p-10,p+10,21)';
    gauss = fittype('a*exp(-((x-m)/s)^2)+c','coeff',{'a','m','s','c'});
    f=fit(y,short,gauss,'StartPoint',[3000, p, 3, me],'Lower',[1000,p-10,1,me*0.5],'Upper',[10000,p+10,4,4*me]);
    coeffs = coeffvalues(f);
    pix(i)=coeffs(2);
    i=i+1;
end
lam = [2.04127;2.05635;2.07291;2.15073;2.18023;2.19556;2.21255;2.23127];

poly3 = fittype('a*x^3+b*x^2+c*x+d','coeff',{'a','b','c','d'});
f=fit(pix,lam,poly3,'StartPoint',[-1.6e-11,8.5e-08,3.4e-4,1.9]);

coeffs = coeffvalues(f);

p0 = coeffs(4);%1.944;
p1 = coeffs(3);%3.035e-4;
p2 = coeffs(2);%1.536e-7;
p3 = coeffs(1);%-5.705e-11;


pixel = linspace(1,1024,1024)';
lambda = p0 + p1*pixel + p2*pixel.^2 + p3*pixel.^3;

%figure;
%plot(lambda,spec);
%save(strcat(path,'lambdaKprime.txt'),'lambda','-ascii');

%lambda1 = dlmread('C:/Users/chrifra/Desktop/BranchingRatio/prepareFits/lookup/lambdaKprime.txt');
end