function [lambda]=ComputeLambda(X,T)

%==========================================================================
% Name: fun_0_4
% Version: 1.1
% Created on: 19/11/2019
% Last modified on: 19/11/2019
%
% Purpose(s): computation of the mixture average thermal conductivity
%
% Inputs:
% - X : array of the species mole fractions
% -	T : temperature
%
% Outputs:
%   - out: computed thermal conductivity
%
% Notes:
%
%==========================================================================
%
%   Copyright (C) 2019 Sorbonne Univ.
%
%==========================================================================

% Global variables
global polyPSC S_PSC mu_PSC SpecNb 

% Mixture average model
lambda=zeros(size(T)); % initialization

X=X+1e-15*ones(SpecNb,1); % to avoid division by 0

for i=1:length(T)

    coeff1=0; coeff2=0;
    
    for k=1:SpecNb

        lambdak=exp(polyval(polyPSC(k,:),log(T(i)),S_PSC(k),mu_PSC(k,:)));
        coeff1=coeff1+X(k)*lambdak;
        coeff2=coeff2+X(k)/lambdak;

    end

    lambda(i)=0.5*(coeff1+1/coeff2);
    
end

    