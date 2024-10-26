%% -----------------Worspace initialization--------------------------------
clc;clear;close all
IsTrace=0;
%% Input
p               = 101325;  % (Pa) ambient pressure
Ti              = 298;     % (K) inlet temperature
grad            = [0.9 0.1];
curv            = [0.9 0.1];
flame_model     = 1;         % 0: burner / 1: free flame model
Energie_solve   = 1;         % 0: temp profile  / 1: energy equation solved
load_former_sol = 0;         % 0: start from scratch / 1: start from former solution
% phi=1;
% H2  = phi;
% O2  = 0.5 ;
% N2  = 0.5*3.76;
% sumM = H2 + O2+ N2 ;
% X_H2 = H2/sumM ;
% X_O2 = O2/sumM ;
% X_N2 = N2/sumM;
% filename = 'resultats1';
% model_opt = [flame_model Energie_solve, load_former_sol];
% REACT_init = [X_H2 X_O2 X_N2];



%% Run
% [debit,Temperature,Yk]=PremixedFlame(REACT_init, p, Ti, grad, curv, model_opt, filename);

%%
% % phi=[0.5 0.61 0.7:0.1:15];
phi=[0.5 0.7 0.9 1.0 1.1 1.3 1.5];
for i = 1 : length(phi)
    H2  = phi(i);
    O2  = 0.5 ;
    N2  = 0.5*3.76;
    sumM = H2 + O2+ N2 ;
    X_H2 = H2/sumM ;
    X_O2 = O2/sumM ;
    X_N2 = N2/sumM;
    REACT_INIT= [X_H2 X_O2 X_N2]
    filename = ['resul_phi' num2str(phi(i)*10,'%02i')]; 
    [debit,Temperature,Yk,GRID]= PremixedFlame(REACT_init, p, Ti, grad, curv, model_opt, filename);
    debit.phi_05
    debit.(['phi_' num2str(phi(i)*10,'%02i')])
end

