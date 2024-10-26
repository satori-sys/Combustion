function out=ThermoMix(meth,X,T,P,Xsto)

%==========================================================================
% Name: ThermoMix
% Version: 1.1
% Created on: 19/11/2019
% Last modified on: 19/11/2019
% Author(s): G. Legros
%
% Purpose(s): computation of the mixture thermodynamic properties
%
% Inputs:
% - meth: name of the property to be computed
% - X: mixture composition [-] (size=size(SpecAtomNam))
% - T: temperature [K] 
% - P: pressure [Pa]
% - Xsto: stoichiometric coefficients (size=size(SpecAtomNam))
%
% Outputs:
%   - out: computed property
%
% Notes:
%
%==========================================================================
%
%   Copyright (C) 2019 Sorbonne Univ.
%
%==========================================================================

%-----------------Global variables ----------------------------------------
global SpecmMolMas
global ThermAhigh ThermAlow
global R P0

posH=find(T>1000);
posL=find(T<=1000);

switch meth
    
    case 'Mm' % Molecular weigh [kg/mole]
        out=X'*SpecmMolMas;
        
    case 'r' % gas constant [J/kg/K]
        out=R/ThermoMix('Mm',X,T,P);

    case 'CpEsp' % specific capacity at constant pressure per unit mol of the pure species [J/mole/K]
        Cpcoef=[1*ones(1,length(T)) ; T ;T.^2 ;T.^3 ;T.^4];
        CpEspeH=R*ThermAhigh(:,1:5)*Cpcoef;
        CpEspeL=R*ThermAlow(:,1:5)*Cpcoef;
        CpEspe=[CpEspeL(:,posL) CpEspeH(:,posH)];
        out=CpEspe;
        
    case 'cpEsp' % specific capacity at constant pressure per unit mass of the pure species [J/kg/K]
        Cpcoef=[1*ones(1,length(T)) ; T ;T.^2 ;T.^3 ;T.^4];
        CpEspeH=R*ThermAhigh(:,1:5)*Cpcoef;
        CpEspeL=R*ThermAlow(:,1:5)*Cpcoef;
        CpEspe=[CpEspeL(:,posL) CpEspeH(:,posH)];
        out=CpEspe./SpecmMolMas;
        
    case 'Cpmo' % specific capacity at constant pressure per unit mole of mixture [J/mole/K]
        Cpcoef=[1*ones(1,length(T)) ; T ;T.^2 ;T.^3 ;T.^4];
        CpEspeH=R*ThermAhigh(:,1:5)*Cpcoef;
        CpEspeL=R*ThermAlow(:,1:5)*Cpcoef;
        CpEspe=[CpEspeL(:,posL) CpEspeH(:,posH)];
        out=X'*CpEspe;
        
    case 'Cpm'  % specific capacity at constant pressure per unit mass of mixture [J/kg/K]
        out=ThermoMix('Cpmo',X,T,P)/ThermoMix('Mm',X,T,P);
        
    case 'Cvmo' % specific capacity at constant volume per unit mole of mixture [J/mole/K]
        Cpmo=ThermoMix('Cpmo',X,T,P);
        out=Cpmo-R;
        
    case 'Cvm' % specific capacity at constant volume per unit mass of mixture [J/kg/K]
        Cpm=ThermoMix('Cpm',X,T,P);
        out=Cpm-ThermoMix('r',X,T,P);
        
    case 'g' % gamma ratio of specific capacities 
        Cpm=ThermoMix('Cpmo',X,T,P);
        Cvm=ThermoMix('Cvmo',X,T,P);
        out=Cpm./Cvm;
        
    case 'htmo' % total enthalpy per unit mole of mixture [J/mole]
        hcoef=[T ;T.^2/2 ;T.^3/3 ;T.^4/4;T.^5/5;1*ones(1,length(T))];
        hEspeH=R*ThermAhigh(:,1:6)*hcoef;
        hEspeL=R*ThermAlow(:,1:6)*hcoef;
        hEspe=[hEspeL(:,posL) hEspeH(:,posH)];
        out=X'*hEspe;

    case 'htm' % total enthalpy per unit mass of mixture [J/kg]
        out=ThermoMix('htmo',X,T,P)/ThermoMix('Mm',X,T,P);
    
    case 'hfmo' % enthalpy of formation per unit mole of mixture [J/mole]
        out=ThermoMix('htmo',X,298,P);
        
    case 'hfm' % enthalpy of formation per unit mass of mixture [J/kg]
        out=ThermoMix('htmo',X,298,P)/ThermoMix('Mm',X,T,P);
        
    case 'hfmT' % sensible enthalpy per unit mole of each species [J/mol]
        hcoef=[T ;T.^2/2 ;T.^3/3 ;T.^4/4;T.^5/5;1*ones(1,length(T))];
        hEspeL=(R*ThermAlow(:,1:6)*hcoef);
        hEspeH=(R*ThermAhigh(:,1:6)*hcoef);
        hEspe=[hEspeL(:,posL) hEspeH(:,posH)];
        out=hEspe;
        
    otherwise
        disp('Unknown method.')
        
end
