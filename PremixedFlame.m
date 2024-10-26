function [debit,Temperature,Yk,GRID]=PremixedFlame(REACT_init, p, T_init, grad, curv, model_opt, filename)


% Created on: 19/11/2019
% Last modified on: 19/11/2019
%
% Purpose(s): main program to resolve the structure of the 1D premixed
% flame
%
% Inputs:
% -
%
% Outputs:
% -
%
% Notes:
%
%==========================================================================
%
%   Copyright (C) 2019 Sorbonne Univ.
%
%==========================================================================

%% -----------------Worspace initialization--------------------------------
% clc
% clear
% clear global
% close all
warning('off'); % avoid warnings display
%%

%% -----------------Global variables definition----------------------------
GlobalData
%%

%% -----------------Parameters definition----------------------------------
param
%%

%--------------------------------------------------------------------------
% Input parameters
%--------------------------------------------------------------------------

%-----------------Options for visualizing the computation progress---------
IsTrace=0; % for visualizing on the prompt
TrackUnst=1; % for writing the unsteady solutions into the output file

%---------Option for loading the parameters of a former computation--------
IsReStart=0; % must the computation start from a former solution
%iSol=1; % 1 for uploading the coarsest former solution
        % 2 for uploading the finest former solution

% if IsReStart
%     
%     FileReStart='out2.txt'; %name of the file containing the former solution
%     FileOut='out.txt'; %name of the output file
%     
%     [free,NRJ,wind,P0,PRES,Ti,FLRT,...
%     nFileSchemCK,nFileThermData,nFileTranspData,...
%     mLow,mHigh,TLow,THigh,YLow,YHigh,MaxItNb,MaxIter,MaxSU,JacAge,JacAgeUnst,...
%     tol_ss,tol_ts,timestep,nbTimeSteps,time,TDEC,TINC,freqinc,dtmax,...
%     MaxTimeStepChange,npts,XSTR,XEND,XCEN,WMIX,TFIX,raffine,ada,...
%     gradmax,curvmax,mod_trans,addVcor,...
%     thermdif,REAC]=ReadHeadReStart(FileReStart);
%     
% else
   
%     [free,NRJ,wind,P0,PRES,Ti,FLRT,FileOut,...
%     nFileSchemCK,nFileThermData,nFileTranspData,nFileTempInit...
%     mLow,mHigh,TLow,THigh,YLow,YHigh,MaxItNb,MaxIter,MaxSU,JacAge,JacAgeUnst...
%     tol_ss,tol_ts,timestep,nbTimeSteps,time,TDEC,TINC,freqinc,dtmax,...
%     MaxTimeStepChange,npts,XSTR,XEND,XCEN,WMIX,TFIX,...
%     raffine,ada,gradmax,curvmax,mod_trans,addVcor,thermdif]=InitComp;

[FLRT,FileOut,...
nFileSchemCK,nFileThermData,nFileTranspData,nFileTempInit...
MaxItNb,MaxIter,MaxSU,JacAge,JacAgeUnst...
timestep,nbTimeSteps,time,freqinc,dtmax,...
MaxTimeStepChange,XSTR,XEND,XCEN,WMIX,...
raffine,ada,gradmax,curvmax]=InitComp;

%--------------------------------------------------------------------------
% Load operator options (2)
%--------------------------------------------------------------------------

gradmax     = grad;
gradmax     = curv;
FileOut     = filename;
free        = model_opt(1); % 0: burner        / 1: free flame model
NRJ         = model_opt(2); % 0: temp profile  / 1: energy equation solved
PRES        = p;            % (Pa) ambient pressure
Ti          = T_init;       % (K)  inlet temperature
TFIX        = T_init+100;          % [K] temperature at x=XCEN for the 1st iteration
    
%--------------------------------------------------------------------------
% Workspace definition
%--------------------------------------------------------------------------
load('FlameVar');

% repert=change_dir;
% 
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % Input files interpretation
% %--------------------------------------------------------------------------
% 
% %-----------------Chemical kinetics data-----------------------------------
% nomfileCK=fullfile(repert.dir,repert.dataSch,nFileSchemCK);
% [SpecAtomNam,SpecNb,ReactNb,ReactCoefStoe,ReactTyp,ReactAbE,ReactTroe,ReactCoefEsp,IsDup,err]=interpkin(nomfileCK);
% 
% %-----------------Thermodynamics data--------------------------------------
% nomfileNASA=fullfile(repert.dir,repert.dataSch,nFileThermData);
% [SpecmMolMas,ThermAhigh,ThermAlow,OriginThermoData,err1,err2]=interpthermo(nomfileCK,nomfileNASA,SpecAtomNam);
% 
% %-----------------Transport data-------------------------------------------
% nomfileTran=fullfile(repert.dir,repert.dataTrans,nFileTranspData);
% [TranCoeff, OriginTranData, err3]= interptransp (SpecAtomNam,nomfileTran);
% [TranCoeff, OriginTranData]= assigntransp (SpecmMolMas,TranCoeff,OriginTranData);
% 
% %-----------------Collision integrals data--------------------------------- 
% nFileCollData1='11.txt';
% nFileCollData2='22.txt';
% nFileCollDataA='A.txt';
% nFileCollDataB='B.txt';
% nFileCollDataC='C.txt';
% nomfileColl=fullfile(repert.dir,repert.dataTrans,nFileCollData1);
% [delta1, T1, collision11, err1]= readcollision (nomfileColl);
% nomfileColl=fullfile(repert.dir,repert.dataTrans,nFileCollData2);
% [delta2, T2, collision22, err2]= readcollision (nomfileColl);
% nomfileColl=fullfile(repert.dir,repert.dataTrans,nFileCollDataA);
% [deltaA, TA, collisionA, errA]= readcollision (nomfileColl);
% nomfileColl=fullfile(repert.dir,repert.dataTrans,nFileCollDataB);
% [deltaB, TB, collisionB, errB]= readcollision (nomfileColl);
% nomfileColl=fullfile(repert.dir,repert.dataTrans,nFileCollDataC);
% [deltaC, TC, collisionC, errC]= readcollision (nomfileColl);
% ind=find(abs([err1 err2 errA errB errC]-ones(1,5))==1);
% ind=[];chaine=[];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Polynomial interpolation for the transport properties
%--------------------------------------------------------------------------

%-----------------Pure species polynomial coefficients---------------------
% lowerT:  [K] lower temperature where the polynoms are defined
% upperT:  [K] upper temperature where the polynoms are defined
% polyPSV(k,:) : polynomial coefficients for the kth species shear vicosity computation
% S_PSV(k,:) and mu_PSV(k:,:): parameters for the evluation of it
% errPSV: relative error at each T due to the interpolation
% polyBDC(k,:), S_BDC(k,:) et mu_BDC(k:,:): same for the binary diffusion
% polyPSC(k,:), S_PSC(k,:) et mu_PSC(k:,:): same for the thermal conductivities
% lowerT=300;
% upperT=5000;
% switch mod_trans
%     case 'MIXTAV' % MIXTURE AVERAGE option
%  
%         [polyPSV,S_PSV,mu_PSV,errPSV,polyBDC,S_BDC,mu_BDC,errBDC,polyPSC,S_PSC,mu_PSC,errPSC,polyTDC,S_TDC,mu_TDC,errTDC]=...
%             polytran(TranCoeff,lowerT,upperT,delta1,T1,collision11,delta2,T2,collision22,deltaA,TA,collisionA,deltaB,TB,collisionB,deltaC,TC,collisionC);
%                     
% %   case 'MULTICOMP' % MULTICOMPONENT option
% %         
%     otherwise
% 
%         disp('Unknown model of transport properties computations...')
%         
% end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%-----------Writing the functions depending on the input parameters--------
%message1=SchemeWrite2(); % reaction rates
InitReStart; % output writing

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
% if IsReStart
% 
%     [x0]=ReadReStart(FileReStart,iSol);
%     
%     UNST=0;       % is the computation unsteady
%     SuccesUnst=0; % was the latest unsteady computation successfull
%     dt=timestep;  % time step initialization
% 
%     LastInc=0;    % increasing (1) or decreasing (0) time step
% 
%     itnb=0;       % number of switch steady/unsteady    
%     nbdt=0;       % number of unsteady computations performed
% 
%     nbref=1;      % number of refinement procedures
% 
% 
% else
%     
    [x0,REAC]=InitSol(FLRT,XSTR,XEND,XCEN,WMIX,REACT_init,T_init,model_opt);
    
    UNST=0;       % is the computation unsteady
    SuccesUnst=0; % was the latest unsteady computation successfull
    dt=timestep;  % time step initialization

    LastInc=0;    % increasing (1) or decreasing (0) time step
    
    itnb=0;       % number of switch steady/unsteady
    nbdt=0;       % number of unsteady computations performed
    
    nbref=1;      % number of refinement procedures
    
%end

%--------------------------------------------------------------------------
% Insert here the parameters to be changed once a former solution uploaded
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%-----------------Head of the output file----------------------------------
[FID1,message]=fopen(FileOut,'wt+');
HeadOut(IsReStart,nFileSchemCK,nFileThermData,nFileTranspData,...
    MaxItNb,MaxIter,MaxSU,JacAge,JacAgeUnst,timestep,...
    nbTimeSteps,time,freqinc,dtmax,MaxTimeStepChange,...
    XSTR,XEND,XCEN,WMIX,raffine,ada,gradmax,curvmax,FLRT,REAC);
fclose(FID1);    

%--------------------------------------------------------------------------
% Iterations
%--------------------------------------------------------------------------

GoOut2=0; %ask for going out from the adaptation loop

[FID1,message]=fopen(FileOut,'at+'); % write into the output file
WriteReStart(x0,time,itnb,dt);

while ~GoOut2
    
    GoOut1=0; %ask for going out from the steady solution search loop

    %-----------------Vectors of bounds creation---------------------------
    VecBound(mLow,mHigh,TLow,THigh,YLow,YHigh);

    while (~GoOut1)&&(itnb<MaxItNb)

        itnb=itnb+1;

        %-----------------Current loop parameter---------------------------
        if ((itnb==1)&(~IsReStart))
            burner=1;
            energie=0; %not solving the energy equation (given temperature profile)
        else
            burner=~free;
            energie=NRJ;
            nbdt=0; % initialization of the number of unsteady computations
        end

        %-----------------Boundary conditions on the left (fresh)----------
        if burner==1
            LBC(1)=FLRT; % in term of flow rate in any case for the 1st iteration
        else
            LBC(1)=x0(1);
        end
        LBC(2)=Ti; % in temperature
        LBC(3:2+SpecNb)=mass2mol('mass',REAC); % in mass fraction

        %-----------------Steady/Unsteady computations---------------------
        itSU=0;
        IsNextStat=1;
        SuccesSteady=0; % was the latest steady computation successfull

        while ((itSU<MaxSU)&(~SuccesSteady))

            itSU=itSU+1; % number of Steady/Unsteady switch 

            %----- solve the steady model provided by the finite differences
            %----- scheme
            if IsNextStat
                %fprintf('Beginning of the steady solution search...\n')
                if IsTrace
                    fprintf(' it.   function norm       step norm     trial value (eqnb)      trial step (vanb)    damp. coeff.   Jac Cond numb\n')
                end
                UNST=0;
                [x,fx,SuccesSteady] = dampednewton('eq_conserv',x0,MaxIter,JacAge);
            end

            %----- solve the unsteady model provided by the finite differences
            %----- scheme
            if ~SuccesSteady

                %save the initial solution
                xinit=x0;

                fx0 = feval('eq_conserv',x,xinit);
                nfx0=log10(norm(fx0,inf));                
                
                %fprintf('Beginning of the unsteady solution search...\n')
                
                if IsTrace
                    fprintf('step       time step   function norm     sol. change   Jac Cond numb\n')
                    fprintf('Steady function norm: %e\n', nfx0)
                end

                UNST=1;
                j=0;
                GoOut1=0;
                Jak=[];
                Age=0;
                LastInc=0;
                kmax=MaxTimeStepChange;

                while ((~GoOut1)&(j<nbTimeSteps))

                    j=j+1;

                    nbdt=nbdt+1;

                    %time step
                    if ((mod((nbdt-1),freqinc)==0)&&((nbdt-1)~=0))
                        dt=min(dtmax,dt*TINC);
                    end

                    GoOut2=0;
                    k=0;

                    while ((~GoOut2)&(k<kmax))

                        [x,fx,SuccesUnSteady,Jak,Age] = dampednewtonUNST('eq_conserv',x0,dt,MaxIter,JacAgeUnst,Jak,Age);

                        UNST=0;
                        fx0 = feval('eq_conserv',x);
                        UNST=1;

                        if IsTrace
                            fprintf('%4d  %14.6e  %14.6e  %14.6e  %14.6e\n',...
                                nbdt, log10(dt), log10(norm(fx0,inf)), log10(norm(x-x0,inf)), log10(cond(Jak)))
                        end
                        
                        if SuccesUnSteady
                            GoOut2=1;
                            Age=Age+1;
                            time=time+dt;
                        else
                            k=k+1;
                            if LastInc
                                dt=dt*TINC;
                            else
                                dt=dt/TDEC;
                            end
                            Jak=[];
                            Age=0;
                            if IsTrace
                                fprintf('new time step: %e \n',log10(dt))
                            end
                        end
                    end

                    if GoOut2
                        x0=x;
                    else
                        if IsTrace
                            fprintf('divergence at time step %d\n',j)
                        end
                    end

                    if ((k==kmax)&&(~LastInc))
                        LastInc=1;
                        kmax=2*MaxTimeStepChange;
                        if IsTrace
                            fprintf('minimum time step => time step increase \n')
                        end
                        j=j-1;
                        nbdt=nbdt-1;
                    end

                    if ((k==kmax)&&(LastInc))
                        if IsTrace
                            fprintf('ultimate divergence at time step %d\n',j)
                        end
                        GoOut1=1;
                        Jak=[];
                    end                    

                end

                %check for the solution evolution along the time stepping
                dx0=x0-xinit;
                evol=sum(abs(dx0));
                if evol~=0

                    %----- regularization before incorporation into the
                    %----- following iteration
                    x0=x';

                    if size(x0,2)>size(x0,1)
                        x0=x0';
                    end

                    if j==nbTimeSteps
                        IsNextStat=1;
                        if IsTrace
                            fprintf('Time stepping fully successfull until time step %d \n', nbdt)
                        end
                        Jak=[];
                    else
                        IsNextStat=1;
                        if IsTrace
                            fprintf('Time stepping partially successfull until time step %d \n', nbdt)
                        end
                        Jak=[];
                    end
                    
                    if TrackUnst
                        WriteReStart(x0,time,itnb,dt);
                    end
                    
                else
                    IsNextStat=0;
                    if IsTrace
                        fprintf('Time stepping failed\n')
                    end
                    Jak=[];
                    GoOut1=1;
                    fprintf('END WITHOUT CONVERGENCE\n');
                    fprintf(FID1,'END WITHOUT CONVERGENCE\n');
                    fprintf(FID1,'END OF FILE\n');
                    fclose(FID1);

                end
            end

        end

        if SuccesSteady
            s1=''; s2='';
            if burner
                s1='BURNER FLAME';
            else
                s1='FREELY PROPAGATING FLAME';
                GoOut1=1; %ask for going out
                itnb=1;
            end
            if energie
                s2='with energy equation solved';
            else
                s2='with fixed temperatures';
            end
            fprintf('Convergence at stage %s %s \n',s1,s2)
            %----- regularization before incorporation into the
            %----- following iteration
            x0=x';

            if size(x0,2)>size(x0,1)
                x0=x0';
            end

            WriteReStart(x0,time,itnb,dt);

        end

    end

    if raffine
        
        GoOut3=0;
        
        while ~GoOut3

            gradmax0=gradmax(nbref);
            curvmax0=curvmax(nbref);
            
            % look for the locations of new mesh points
            ind_ref=mustref(GRID,x0,ada,gradmax0,curvmax0);

            if ~isempty(ind_ref)
                
                GoOut2=0;

                [debit0,Temperature0,Yk0]=de_composition(x0);

                %new mesh points location
                newpt=[]; ord=[]; GRID0=[];
                GRID0=GRID;
                newpt=(GRID0(ind_ref)+GRID0(ind_ref+1))/2;
                [GRID,ord]=sort( [GRID0; newpt]);

                %interpolation for the variables values at the new mesh points
                newdebit = [];
                newdebit = interp1(GRID0,debit0,newpt,'linear');
                debit0 = [debit0 newdebit'];
                debit0 = debit0(ord);

                newTemp = [];
                newTemp = interp1(GRID0,Temperature0,newpt,'linear');
                Temperature0 = [Temperature0 newTemp'];
                Temperature0 = Temperature0(ord);

                newY = [];
                newY = interp1(GRID0,Yk0',newpt,'linear');
                Yk0 = [Yk0 newY'];
                Yk0(:,:) = Yk0(:,ord);

                %new index of TFIX
                jfix=find(Temperature0==TFIX);

                [x0]=re_composition(debit0,Temperature0,Yk0);

                %time step initialization
                dt=timestep;
                
                GoOut3=1;

            else
                if nbref==length(gradmax) %when the most severe conditions is fullfilled
                    GoOut3=1;
                    GoOut2=1;
                    fprintf('SUCCESSFULL END\n');
                    fprintf(FID1,'SUCCESSFULL END\n');
                    fprintf(FID1,'END OF FILE\n');
                    fclose(FID1);
                else
                    nbref=nbref+1;
                end
            end
            
        end
        
    else
        
        GoOut2=1;
        
    end

end

%--------------------------------------------------------------------------
% Affichage
%--------------------------------------------------------------------------
[debit,Temperature,Yk]=de_composition(x);
% figure
% plot(GRID,debit)
% xlabel('x (cm)'); ylabel('debit (g/cm²/s)');
% title('debit')
% figure
% plot(GRID,Temperature)
% xlabel('x (cm)'); ylabel('T (K)');
% title('temperature')
% figure
% plot(GRID,Yk(3,:),'b',GRID,Yk(9,:),'g',GRID,Yk(4,:),'r',GRID,100*Yk(6,:),'k')
% legend('H2','O2','H2O','100 x HO2')
% xlabel('x (cm)'); ylabel('Y_k (-)');



end