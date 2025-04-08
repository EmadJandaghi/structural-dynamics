function [DM,VM,AM,aTime]=MDEDNL(gamma,M,C,activeDOFs,time,dt,D0,V0)
%% Info
% Based on:
% has not been published yet


% a Model dependent time integration technique with effective numerical
% damping. (MDED)
% Based on the implicit bathe technique

% Written by: A. H. Namadchi (1/2018)
%% Core
problemData;
% global P
P0=P(0);
nDOFs=length(P0);
% [K,fint0]=MindlinPlateKtFint(D0,elConn,DOFConn,...
%     refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
[K,fint0]=T3DKtFint(D0,elConn,DOFConn,refCoords,areaList...
    ,L0List,s0List,activeDOFs,E,1);

A0=M\(P0(activeDOFs)-fint0-C*V0);
clear P0

totalSteps=length(0:dt:time);
[DM,VM,AM]=deal(zeros(nDOFs,totalSteps));

DM(activeDOFs,1)=D0;
VM(activeDOFs,1)=V0;
AM(activeDOFs,1)=A0;

% Preliminaries
j=1;
nADOFs=length(activeDOFs);
IM=eye(nADOFs);

% Weighting Factors -------------------------------------------------------
g0=(-2+gamma)^2;g1=2*(-2+gamma)*(-2+gamma^2);
g2=2*(-2+gamma)*(-1+gamma)*gamma;g3=2-2*gamma+gamma^2;
g4=(-1+gamma)*gamma*(-2+gamma^2);g5=(-1+gamma)*gamma;
g6=(-2+gamma)*(-1+gamma);g7=(-1+gamma)^2;g8=4-gamma^2*(6+(-4+gamma)*gamma);
g9=(1-gamma)*gamma*(2-2*gamma+gamma^2);g10=2*(-2+gamma)*(2-2*gamma+gamma^2);
g11=-6+gamma*(4+gamma);g12=-2+gamma^2;
% -------------------------------------------------------------------------
% Integration Matrices
MNAlpha0=4*g0*M + g1*dt*C - g9*(dt^2)*K ...
         + g2*(dt^2)*(C*(M\C)) + (g5^2)*(dt^3)*(C*(M\K));
MNAlpha1=2*g0*M + g2*dt*C + (g5^2)*(dt^2)*K;  
MTheta0=4*g0*M + g1*dt*C + (g3^2)*(dt^2)*K +...
        g2*(dt^2)*(C*(M\C)) + g4*(dt^3)*(C*(M\K)) + (g5^2)*(dt^4)*(K*(M\K));
    
alpha3=(4*g0*IM + g2*dt*(M\C) - g9*(dt^2)*(M\K))\(2*g0*IM);
alpha4=(4*g0*IM + g2*dt*(M\C) - g9*(dt^2)*(M\K))\...
       (2*g0*IM + g2*dt*(M\C) + (g5^2)*(dt^2)*(M\K));
inMTheta2=(4*g0*M + g2*dt*(C) - g9*(dt^2)*(K))\IM;
% Factorizations
% I didn't use inv(A) since it's less accurate
inKeq=MTheta0\IM;
inMeq=(M+dt*C*alpha4)\IM;
kInvM=K/M;
cInvM=C/M;
tic;
for t=dt:dt:time
    P1=P(t);
    Pg=P(t-gamma*dt);
    P0=P(t-dt);
    NAlpha2=( (g11-4*g5)*IM - (g2+g5)*dt*cInvM - (dt^2)*(g5^2)*kInvM )*P0(activeDOFs)+...
            ( -g12*IM - g5*dt*cInvM )*Pg(activeDOFs)+...
            ( 4*g7*IM + (g2+2*g5)*dt*cInvM + (g5^2)*(dt^2)*kInvM )*P1(activeDOFs);
    NAlpha5=( -2*g6*IM )*P0(activeDOFs)+...
            ( (g3-g12)*IM )*Pg(activeDOFs)+...
            ( (3*g12-g11)*IM )*P1(activeDOFs);
    deltaD=inKeq*(dt*MNAlpha0*V0+(dt^2)*(MNAlpha1)*A0+(dt^2)*NAlpha2);
    D1=D0+deltaD;
    fint=T3DFint(D1,elConn,DOFConn,refCoords,areaList...
                ,L0List,s0List,activeDOFs,E,1);
    A1=inMeq*(P1(activeDOFs)-fint-...
              C*(V0+dt*(alpha3)*A0+dt*(inMTheta2*(NAlpha5))));
    V1=V0+dt*(alpha3)*A0+dt*(alpha4)*A1+dt*(inMTheta2*(NAlpha5));
    
    j=j+1;
    sprintf('MDED NLinear Progress: %0.2f %',(j/totalSteps)*100)
    DM(activeDOFs,j)=D1;
    VM(activeDOFs,j)=V1;
    AM(activeDOFs,j)=A1;
    D0=D1;V0=V1;A0=A1;
    
end
aTime=toc;

end