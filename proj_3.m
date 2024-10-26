% Exercice 3
% workspace preparation
clear all
clc
close all

% load the required variables
load('FlameVar.mat');
% load('CoeffThermo.mat');

phi_05 = readtable('phi_05.xlsx','ReadRowNames',false);
phi_07 = readtable('phi_07.xlsx','ReadRowNames',false);
phi_09 = readtable('phi_09.xlsx','ReadRowNames',false);
phi_10 = readtable('phi_10.xlsx','ReadRowNames',false);
phi_11 = readtable('phi_11.xlsx','ReadRowNames',false);
phi_13 = readtable('phi_13.xlsx','ReadRowNames',false);
phi_15 = readtable('phi_15.xlsx','ReadRowNames',false);

% Récupération des données
X_05 = phi_05.X;
X_07 = phi_07.X;
X_09 = phi_09.X;
X_10 = phi_10.X;
X_11 = phi_11.X;
X_13 = phi_13.X;
X_15 = phi_15.X;

T_05 = phi_05.T;
T_07 = phi_07.T;
T_09 = phi_09.T;
T_10 = phi_10.T;
T_11 = phi_11.T;
T_13 = phi_13.T;
T_15 = phi_15.T;

% Figures

% Question 1
plot(X_05,T_05);
hold on 
plot(X_07,T_07);
plot(X_09,T_09);
plot(X_10,T_10);
plot(X_11,T_11);
plot(X_13,T_13);
plot(X_15,T_15);
hold off
title('Température en fonction de la position pour phi[0.5,0.7,0.9,1.0,1.1,1.3,1.5]')
xlabel('Position (pt)')
ylabel('Température T[K]')
legend('phi=0.5','phi=0.7','phi=0.9','phi=1.0','phi=1.1','phi=1.3','phi=1.5')

% Question 2
H2O = phi_10.H2O;
X = phi_10.X;
V = phi_10.V;

figure;
plot(X,H2O)
xlabel('Position x[m]')
ylabel('Fraction molaire X_{H2O} [-]')
title('Fraction molaire')

figure;
plot(X,V)
xlabel('Position x[m]')
ylabel('Vitesse de réaction V_{H2O} [m.s-1]')
title('Vitesse de réaction')

% Exercice Bonus
img1 = imread('1.png');
img2 = imread('2.png');
img3 = imread('3.png');
img4 = imread('4.png');

myVideo = VideoWriter('myvideo.avi');
myVideo.FrameRate = 30;
open(myVideo);
% % resize the image to 914 by 744
% img_resized = imresize(img2, [744 914]);
% % save the resized image
% imwrite(img_resized, 'img2re.jpg');

writeVideo(myVideo, img1);
writeVideo(myVideo, img2);
writeVideo(myVideo, img3);
writeVideo(myVideo, img4);

close(myVideo);

% Récupération des données
% a1=[0,00000000e+00
% 1,66656608e-01
% 3,00000000e-01
% 6,00000000e-01
% 9,00000000e-01
% 1,20000000e+00
% 1,50000000e+00
% ]
% b1=[2,98000000e+02
% 3,98000000e+02
% 5,34749046e+02
% 1,90090090e+03
% 1,96599825e+03
% 1,97297297e+03
% 1,97297297e+03
% ]
% plot(a1,b1)
%xm1=[0,00000000e+00
% 1,67E-01
% 3,00E-01
% 6,00E-01
% 9,00E-01
% 1,20E+00
% 1,50E+00
% ]
% ym1=[2,98000000e+02
% 3,98E+02
% 5,35E+02
% 1,90E+03
% 1,97E+03
% 1,97E+03
% 1,97E+03
% ]
% xm2=[0,00000000e+00
% 8,33E-02
% 1,25E-01
% 1,46E-01
% 1,56E-01
% 1,67E-01
% 1,75E-01
% 1,83E-01
% 1,92E-01
% 2,00E-01
% 2,17E-01
% 2,33E-01
% 2,50E-01
% 2,67E-01
% 3,00E-01
% 4,50E-01
% 6,00E-01
% 7,50E-01
% 9,00E-01
% 1,20E+00
% 1,50E+00
% ]
% ym2=[2,98000000e+02
% 2,98E+02
% 2,98E+02
% 3,01E+02
% 3,13E+02
% 3,98E+02
% 6,62E+02
% 1,01E+03
% 1,30E+03
% 1,50E+03
% 1,68E+03
% 1,78E+03
% 1,86E+03
% 1,91E+03
% 1,98E+03
% 2,11E+03
% 2,18E+03
% 2,23E+03
% 2,26E+03
% 2,29E+03
% 2,29E+03
% ]
% xm4=[0,00000000e+00
% 8,33E-02
% 1,25E-01
% 1,46E-01
% 1,56E-01
% 1,61E-01
% 1,67E-01
% 1,71E-01
% 1,75E-01
% 1,83E-01
% 1,92E-01
% 2,00E-01
% 2,17E-01
% 2,33E-01
% 2,50E-01
% 2,67E-01
% 3,00E-01
% 4,50E-01
% 6,00E-01
% 7,50E-01
% 9,00E-01
% 1,20E+00
% 1,50E+00
% ]
% ym4=[2,98000000e+02
% 2,98E+02
% 2,98E+02
% 3,00E+02
% 3,08E+02
% 3,28E+02
% 3,98E+02
% 5,35E+02
% 7,27E+02
% 1,11E+03
% 1,38E+03
% 1,55E+03
% 1,71E+03
% 1,81E+03
% 1,88E+03
% 1,93E+03
% 2,00E+03
% 2,12E+03
% 2,18E+03
% 2,23E+03
% 2,25E+03
% 2,29E+03
% 2,29E+03
% ]
% xm3=[0,00000000e+00
% 8,33E-02
% 1,25E-01
% 1,46E-01
% 1,51E-01
% 1,56E-01
% 1,59E-01
% 1,60E-01
% 1,61E-01
% 1,62E-01
% 1,63E-01
% 1,63E-01
% 1,64E-01
% 1,65E-01
% 1,65E-01
% 1,66E-01
% 1,67E-01
% 1,67E-01
% 1,68E-01
% 1,68E-01
% 1,69E-01
% 1,69E-01
% 1,70E-01
% 1,70E-01
% 1,71E-01
% 1,71E-01
% 1,72E-01
% 1,72E-01
% 1,73E-01
% 1,73E-01
% 1,74E-01
% 1,74E-01
% 1,75E-01
% 1,76E-01
% 1,77E-01
% 1,78E-01
% 1,79E-01
% 1,80E-01
% 1,81E-01
% 1,82E-01
% 1,83E-01
% 1,84E-01
% 1,85E-01
% 1,86E-01
% 1,87E-01
% 1,89E-01
% 1,90E-01
% 1,91E-01
% 1,92E-01
% 1,93E-01
% 1,94E-01
% 1,95E-01
% 1,96E-01
% 1,97E-01
% 1,98E-01
% 1,99E-01
% 2,00E-01
% 2,02E-01
% 2,04E-01
% 2,06E-01
% 2,08E-01
% 2,12E-01
% 2,17E-01
% 2,25E-01
% 2,33E-01
% 2,42E-01
% 2,50E-01
% 2,67E-01
% 2,83E-01
% 3,00E-01
% 3,38E-01
% 3,75E-01
% 4,50E-01
% 5,25E-01
% 6,00E-01
% 6,75E-01
% 7,50E-01
% 9,00E-01
% 1,05E+00
% 1,20E+00
% 1,50E+00
% ]
% ym3=[2,98000000e+02
% 2,98E+02
% 2,98E+02
% 2,98E+02
% 2,99E+02
% 3,01E+02
% 3,04E+02
% 3,08E+02
% 3,13E+02
% 3,17E+02
% 3,22E+02
% 3,28E+02
% 3,37E+02
% 3,47E+02
% 3,61E+02
% 3,78E+02
% 3,98E+02
% 4,17E+02
% 4,38E+02
% 4,62E+02
% 4,88E+02
% 5,15E+02
% 5,44E+02
% 5,73E+02
% 6,03E+02
% 6,34E+02
% 6,65E+02
% 6,96E+02
% 7,27E+02
% 7,58E+02
% 7,88E+02
% 8,19E+02
% 8,49E+02
% 9,08E+02
% 9,65E+02
% 1,02E+03
% 1,07E+03
% 1,13E+03
% 1,17E+03
% 1,22E+03
% 1,26E+03
% 1,30E+03
% 1,34E+03
% 1,38E+03
% 1,41E+03
% 1,44E+03
% 1,47E+03
% 1,49E+03
% 1,52E+03
% 1,54E+03
% 1,56E+03
% 1,58E+03
% 1,59E+03
% 1,61E+03
% 1,62E+03
% 1,64E+03
% 1,65E+03
% 1,67E+03
% 1,70E+03
% 1,71E+03
% 1,73E+03
% 1,76E+03
% 1,79E+03
% 1,84E+03
% 1,87E+03
% 1,90E+03
% 1,93E+03
% 1,97E+03
% 2,01E+03
% 2,04E+03
% 2,08E+03
% 2,12E+03
% 2,17E+03
% 2,20E+03
% 2,22E+03
% 2,24E+03
% 2,26E+03
% 2,28E+03
% 2,30E+03
% 2,31E+03
% 2,31E+03
% ]
% 
% plot(xm1,ym1)
% plot(xm2,ym2)
% plot(xm3,ym3)
% 
% 
% % legend({'end of solution_1','end of solution_{x+1<n}','end of solution_{n+1}','end of solution_{n}'},'FontWeight','bold', 'Location', 'BestOutside')
% % xlabel('Position','FontWeight','bold');
% % ylabel('T[K]','FontWeight','bold')
% % title(sprintf('mesh and solution refinement'));
