function [DM,VM,AM,aTime]=BatheNLNR(gamma,M,C,activeDOFs,time,dt,D0,V0,...
    maxIter,tolD)
%% Info
% Based on:
% (2005) - K. J. Bathe,M. I. Baig - On a Composite Implicit Time
% Integration Procedure for Nonlinear Dynamics


% The procedure uses two sub-steps per time step dt. in the first sub-step
% the usual trapezoidal rule is used and in the second sub-step a
% three-point backward difference approximation is used. the algorithm
% remains stable for large time step sizes; but when the time step is too
% large the numerical damping can be appreciable. the composite scheme can
% be significantly more effective than the trapezoidal rule when large
% deformations over long time ranges need be calculated.

% gamma = 0.5 is suggested

% Written by: A. H. Namadchi (9/2017)

%% Core
tic;
problemData;
global P
P0=P(0);
nDOFs=length(P0);

% [~,fint0]=T3DKtFint(D0,elConn,DOFConn,refCoords,areaList...
%     ,L0List,s0List,activeDOFs,E,1);
[~,fint0]=MindlinPlateKtFint(D0,elConn,DOFConn,...
    refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);

A0=M\(P0(activeDOFs)-fint0-C*V0);

clear P0 K0

totalSteps=length(0:dt:time);

[DM,VM,AM]=deal(zeros(nDOFs,totalSteps));

DM(activeDOFs,1)=D0;
VM(activeDOFs,1)=V0;
AM(activeDOFs,1)=A0;

j=1;


for t=dt:dt:time
    P1=P(t);
    Pg=P(t-gamma*dt);
    %% NR for the first substep
    Dg=D0;
%     [Ktg,fintg]=T3DKtFint(Dg,elConn,DOFConn,refCoords,areaList...
%     ,L0List,s0List,activeDOFs,E,1);
    [Ktg,fintg]=MindlinPlateKtFint(Dg,elConn,DOFConn,...
                    refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
	KeqG=(4/((dt*gamma)^2))*M+...
          (2/(dt*gamma))*C+Ktg;
    PeqG=Pg(activeDOFs)-fintg+...
         M*((4/((dt*gamma)^2))*(D0-Dg)+...
            (4/(dt*gamma))*V0+A0)+...
         C*((2/(dt*gamma))*(D0-Dg)+V0);
    % Newton Raphson algorithm
    for n=1:maxIter
        deltaDg=KeqG\PeqG;
        Dg=Dg+deltaDg;
%         [Ktg,fintg]=T3DKtFint(Dg,elConn,DOFConn,refCoords,areaList...
%                     ,L0List,s0List,activeDOFs,E,1);
        [Ktg,fintg]=MindlinPlateKtFint(Dg,elConn,DOFConn,...
                        refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
                    
        KeqG=(4/((dt*gamma)^2))*M+...
              (2/(dt*gamma))*C+Ktg;
        PeqG=Pg(activeDOFs)-fintg+...
             M*((4/((dt*gamma)^2))*(D0-Dg)+...
                (4/(dt*gamma))*V0+A0)+...
             C*((2/(dt*gamma))*(D0-Dg)+V0);
        if norm(deltaDg)<=tolD
            sprintf('NR Converged: %d Iterations',n)
            break
        end        
    end
    %% Post-Convergence         
    Vg=(2/(dt*gamma))*(Dg-D0)-V0;
%     Ag=(4/((dt*gamma)^2))*(Dg-D0)+...
%         (-4/(dt*gamma))*V0-A0;
	%% NR for the second substep
    D1=Dg;
% 	[Kt,fint]=T3DKtFint(D1,elConn,DOFConn,refCoords,areaList...
%         ,L0List,s0List,activeDOFs,E,1);
    [Kt,fint]=MindlinPlateKtFint(D1,elConn,DOFConn,...
                    refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
                
    Keq=(((-2+gamma)/(dt*(-1+gamma)))^2)*M+...
        ((-2+gamma)/(dt*(-1+gamma)))*C+Kt;                    
    Peq=P1(activeDOFs)-fint+...
        M*( ((2-gamma)/(dt^2*(-1+gamma)^2*gamma))*(Dg)-...
            (((-2+gamma)/(dt*(-1+gamma)))^2)*D1+... % Linearized term
            ((-2+gamma)/(dt^2*gamma))*D0+...
            (-(1/(dt*(-1+gamma)*gamma)))*Vg+...
            ((-1+gamma)/(dt*gamma))*V0 )+...
        C*( (-(1/(dt*(-1+gamma)*gamma)))*(Dg)-...
            ((-2+gamma)/(dt*(-1+gamma)))*D1+...     % Linearized term
            ((-1+gamma)/(dt*gamma))*D0);
    % Newton Raphson algorithm
    for n=1:maxIter
        deltaD=Keq\Peq;
        D1=D1+deltaD;
%         [Kt,fint]=T3DKtFint(D1,elConn,DOFConn,refCoords,areaList...
%                         ,L0List,s0List,activeDOFs,E,1);
        [Kt,fint]=MindlinPlateKtFint(D1,elConn,DOFConn,...
                    refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
        Keq=(((-2+gamma)/(dt*(-1+gamma)))^2)*M+...
            ((-2+gamma)/(dt*(-1+gamma)))*C+Kt;                    
        Peq=P1(activeDOFs)-fint+...
            M*( ((2-gamma)/(dt^2*(-1+gamma)^2*gamma))*(Dg)-...
                (((-2+gamma)/(dt*(-1+gamma)))^2)*D1+... % Linearized term
                ((-2+gamma)/(dt^2*gamma))*D0+...
                (-(1/(dt*(-1+gamma)*gamma)))*Vg+...
                ((-1+gamma)/(dt*gamma))*V0 )+...
            C*( (-(1/(dt*(-1+gamma)*gamma)))*(Dg)-...
                ((-2+gamma)/(dt*(-1+gamma)))*D1+...     % Linearized term
                ((-1+gamma)/(dt*gamma))*D0);
        if norm(deltaD)<=tolD
            sprintf('NR Converged: %d Iterations',n)
            break
        end        
    end
    %% Post-Convergence         
    V1=((-2+gamma)/(dt*(-1+gamma)))*D1+...
        (1/(dt*(-1+gamma)*gamma))*Dg+...
        ((1-gamma)/(dt*gamma))*D0;
    A1=((-2+gamma)^2/(dt^2*(-1+gamma)^2))*D1+...
        ((-2+gamma)/(dt^2*(-1+gamma)^2*gamma))*Dg+...
        ((2-gamma)/(dt^2*gamma))*D0+...
        (1/(dt*(-1+gamma)*gamma))*Vg+...
        ((1-gamma)/(dt*gamma))*V0;
    
    
    j=j+1;
    sprintf('Bathe Nonlinear Progress: %0.2f %',(j/totalSteps)*100)
    DM(activeDOFs,j)=D1;
    VM(activeDOFs,j)=V1;
    AM(activeDOFs,j)=A1;
    D0=D1;V0=V1;A0=A1;
    
end
aTime=toc;
end