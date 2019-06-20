close all;
clearvars;

v1i = [9 8 8 7 6 5 4 3];
v2i = [7 6 5 4 4 3 2 1];
T = 200;

PHYS_kB = 1/0.6950356; % in cm^-1/K
B = 20;

M = dlmread('13_hit122018.csv');
nu = M(:,1);
lambda=1./nu*10000;
A = M(:,3);
El = M(:,4);
Eu = El + nu;
gu1 = M(:,5);
gu2 = M(:,6);
v1 = M(:,7);
v2 = M(:,8);
j1 = (gu1/2-1)/2;
j2 = (gu2/2-1)/2;

Q = PHYS_kB.*T/(B);

linecolours = [[0 0 1]; [1 0 0]; [0 0 1]; [1 0 0]; [0 0 1]; [1 0 0]; [0 0 1]; [1 0 0]];

fig = figure;
% subplot('Position',[0.13 0.11 0.84 0.87]);
set(fig,'defaultAxesColorOrder',[[0 0 0];[0 0 0]]);
yyaxis left;
for i = 1 : 8
    lines = (v1==v1i(i)).*(v2==v2i(i)).*((j1==j2-1)+(j1==j2)+(j1==j2+1));%.*(j1<6);
    Ap = A(lines==1);
    Eup = Eu(lines==1);
    Elp = El(lines==1);
    lambdap = lambda(lines==1);
    j1p = j1(lines==1);
    gu1p = gu1(lines==1);
    nup = nu(lines==1);
    
    Ip = (gu1p.*exp(-1.*(Eup-min(Eup))/(PHYS_kB.*T)).*Ap/Q)/70;
    for j = 1 : size(Ip)
        plot([lambdap(j) lambdap(j)],[0 Ip(j)],'-', 'Color', linecolours(i,:)); hold on;
    end
end
axis([1.14 1.42 0 1]);
ylabel('Relative population in upper $v$ state [a.u.]','Interpreter','latex');

yyaxis right;
J = dlmread('NOTCamFilters/1.dat');
plot(J(:,1),J(:,2),'-'); hold on;

H = dlmread('NOTCamFilters/3.dat');
plot(H(:,1),H(:,2),'-');

%Kp = dlmread('NOTCamFilters/5.dat');
%plot(Kp(:,1),Kp(:,2),'--');

K = dlmread('NOTCamFilters/8.dat');
plot(K(:,1),K(:,2),'-');

SABER_K = dlmread('SABER_K_filter.txt');
plot(10000./SABER_K(:,1), 100*SABER_K(:,2), '-');

SABER_H = sortrows(dlmread('SABER_H_band.txt'));
plot(10000./SABER_H(:,1), 100*SABER_H(:,2), '-');


xlabel('Wavelength [$\mu$m]', 'Interpreter', 'latex');
ylabel('Filter Transmission [\%]', 'interpreter', 'latex');
ylim([0,100]);

text(1.17,25,'(7,4)','Color','r');
text(1.24,35,'(8,5)','Color', 'b');
% text(1.48,25,'$(3,1)$','Color','r', 'Interpreter', 'latex');
% text(1.55,40,'$(4,2)$','Color','b', 'Interpreter', 'latex');
% text(1.65,55,'$(5,3)$','Color','r', 'Interpreter', 'latex');
% text(1.78,72,'$(6,4)$','Color','b', 'Interpreter', 'latex');
% text(1.95,77,'$(8,6)$','Color','r', 'Interpreter', 'latex');
% text(2.1,77,'$(9,7)$','Color', 'b', 'Interpreter', 'latex');


text(1.32,85,'J','Color','k');
text(1.75,85,'NOTCam: $J$','Color','k', 'Interpreter', 'latex');
% text(1.75,85,'SABER: $1.6$ $\mu$m','Color','k', 'Interpreter', 'latex');
% text(1.9,85,'K^\prime','Color','k');
% text(2.36,85,'NOTCam: $K$','Color','k', 'Interpreter', 'latex');
% text(1.75,85,'SABER: $2.0$ $\mu$m','Color','k', 'Interpreter', 'latex');


%set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 21 12]);
%print(gcf,'3FilterOverview','-dpng',['-r','1000'],'-opengl')