% workspace preparation
clear;clc
close all

% global variables to be used
% polyPSC: polynomial coefficients for pure species conductivity evaluation
% global polyPSC
% S_PSC: error estimates for pure species conductivity evaluation
% global S_PSC
% mu_PSC: mean and standard deviation for pure species conductivity evaluation
% global mu_PSC

% load the required variables
load FlameVar.mat

% parameter of the study
T=(200:100:3000);
X1=zeros(SpecNb);
X2=zeros(SpecNb);
X3=zeros(SpecNb);
X4=zeros(SpecNb);
X5=zeros(SpecNb);
X6=zeros(SpecNb);
X7=zeros(SpecNb);
X8=zeros(SpecNb);
X9=zeros(SpecNb);
X10=zeros(SpecNb);
X11=zeros(SpecNb);
X12=zeros(SpecNb);

% identification of species among SpecAtomNam
Ind_Ar=1;
Ind_N2=7;
Ind_H2O=4;
Ind_O2=9;

% building X1 (argon)
X1(1)=1;
% building X2 (air)
X2(7)=0.79;
X2(9)=0.21;
X2(4)=0;
% building X3
X3(1)=1;
% building X4
X4(2)=1;
% building X5
X5(3)=1;
% building X6
X6(4)=1;
% building X7
X7(5)=1;
% building X8
X8(6)=1;
% building X9
X9(7)=1;
% building X10
X10(8)=1;
% building X11
X11(9)=1;
% building X12
X12(10)=1;

% computations of lambda as a function of T
lambda1 = ComputeLambda(X1,T)*1e3;
lambda2 = ComputeLambda(X2,T)*1e3;
lambda3 = ComputeLambda(X3,T)*1e3;
lambda4 = ComputeLambda(X4,T)*1e3;
lambda5 = ComputeLambda(X5,T)*1e3;
lambda6 = ComputeLambda(X6,T)*1e3;
lambda7 = ComputeLambda(X7,T)*1e3;
lambda8 = ComputeLambda(X8,T)*1e3;
lambda9 = ComputeLambda(X9,T)*1e3;
lambda10 = ComputeLambda(X10,T)*1e3;
lambda11 = ComputeLambda(X11,T)*1e3;
lambda12 = ComputeLambda(X12,T)*1e3;

% draw the figure 1
figure;
subplot(1,2,1);
plot(T,lambda1,T,lambda2,'LineWidth',2)
legend({'X_{Ar}=1 ','X_{N2}=0.79 / X_{O2}=0.21 / X_{H2O}=0'},'FontWeight','bold', 'Location', 'southoutside');
xlabel('Temperature [K]','FontWeight','bold');
ylabel('thermal conductivity [mW/m/K]','FontWeight','bold');
ylim([0,1000]);
xlim([200,3000]);
title(sprintf('thermal conductivity \nof mixtures as a function of T'));
set(gcf, 'Position', [100 100 1200 600]); % [left, bottom, width, height]
set(gca, 'FontSize', 16); % set font size for axis labels and tick labels
set(gca, 'FontWeight', 'bold'); % set font weight for axis tick labels

% draw the figure 2
subplot(1,2,2);
% plot(T,lambda3,T,lambda4,T,lambda5,T,lambda6,T,lambda7,T,lambda8,T,lambda9,T,lambda10,T,lambda11,T,lambda12,'LineWidth',2)
plot(T, lambda3, 'r', 'LineWidth', 2);
hold on;
plot(T, lambda4, 'g', 'LineWidth', 2);
plot(T, lambda5, 'b', 'LineWidth', 2);
plot(T, lambda6, 'm', 'LineWidth', 2);
plot(T, lambda7, 'c', 'LineWidth', 2);
plot(T, lambda8, 'y', 'LineWidth', 2);
plot(T, lambda9, 'k', 'LineWidth', 2);
plot(T, lambda10, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
plot(T, lambda11, 'Color', [0.7 0.2 0.5], 'LineWidth', 2);
plot(T, lambda12, 'Color', [0.4 0.8 0.1], 'LineWidth', 2);
legend({'Argon','Hydrogen','Dihydrogen',"H_2O","H_2O_2","HO_2","Dinitrogen","Oxygen","Dioxygen","OH"},'FontWeight','bold', 'Location', 'BestOutside')
xlabel('Temperature [K]','FontWeight','bold');
ylabel('thermal conductivity [mW/m/K]','FontWeight','bold');
ylim([0,1000]);
xlim([200,3000]);
title(sprintf('thermal conductivity\n of pure species as a function of T'));
set(gcf, 'Position', [50 100 1500 600]); % [left, bottom, width, height]
set(gca, 'FontSize', 16); % set font size for axis labels and tick labels
set(gca, 'FontWeight', 'bold'); % set font weight for axis tick labels