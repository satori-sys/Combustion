% Projet Combustion
% Exercice 1
% workspace preparation
clear
clear global
close all
clc

% load the required variables
load('FlameVar.mat')

% parameters of the study
global SpecmMolMas
global ThermAhigh ThermAlow
global R P0

% Variables
% Temperature vector [K]
T=[300:100:2500];
% Pression atmosphérique [Pa]
P = 1e5;
% building X fraction molaire [-]
X = [0:0.01:1];

% Argument de sortie doit être affecté de l'enthalpie par unitité de masse du mélange [J/kg] 
meth = 'htm'; 

% identification of species among SpecAtomNam
ind_H2O = 4;
ind_N2 = 7;
ind_O2 = 9;

% Fabrication du mélange
X_mel = zeros(SpecNb,length(X));

for k = 1:length(X_mel)
    X_mel(ind_H2O,k)=X(k);
    X_mel(ind_O2,k)=0.21*(1-X(k));
    X_mel(ind_N2,k)=0.79*(1-X(k));
end

% Matrice Tab1 contenant les valeurs des propriétés thermodynamiques pour chaque valeur de X à différentes températures
Tab1=zeros(size(X_mel,2),size(T,2));

for k = 1:size(X_mel,2)
    Tab1(k,:)=ThermoMix(meth,X_mel(:,k),T,P);
end

% Calcul de la valeur exacte de l'enthalpie totale du mélange X_H2O=0.5, X_air=0.5 et T=1400 K
% X = 0.5 en colonne 6
Val = ThermoMix(meth,X_mel(:,6),1400,P); 

% Figure 
levels=(-15:1:15);
contourf((1-X),T,Tab1.'*1e-6,levels,'ShowText','on','lineColor', 'w');

xlabel('X_{AiR} [-]');
ylabel('T [K]');
title({'Enthalpie massique totale du mélange en fonction de la fraction molaire (MJ/kg)'});

%h1 = axes();
%h1.XDir="reverse";
%colorbar;
h2.LevelList=(-14:3);
colormap('summer')
h2.EdgeColor='white';