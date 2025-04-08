function [DM,VM,AM]=HHTAlphaL(alpha,M,K,C,activeDOFs,time,dt,D0,V0)
%% Info
% Based on:
% (1977) - H. M. Hilber, J.R. Hughes, R. L. Taylor - Improved Numerical
% Dissipation For Time Integration Algorithms in Structural Dynamics

% The Hilber-Hughes-Taylor-? algorithm (HHT-?) has been developed with the
% motivation, that the resulting algorithm is unconditionally stable,
% possesses high frequency numerical dissipation that can be controlled by
% a free parameter rather than by the size of the time step such that it
% does not affect the lower modes too strongly.

% The algorithm is second-order time accurate and
% unconditionally stable for:
% -1/3<=alpha<=0, gamma =0.5-alpha and beta=0.25(1-alpha)^2

%% My Note:
% This version is obtained by setting alphaM=0 in G-alpha method.

% Written by: A. H. Namadchi (9/2017)
%% Core
P0=P(0);
A0=M\(P0(activeDOFs)-K*D0-C*V0);
clear P0

totalSteps=length(0:dt:time);
nDOFs=length(M);
[DM,VM,AM]=deal(zeros(nDOFs,totalSteps));

DM(activeDOFs,1)=D0;
VM(activeDOFs,1)=V0;
AM(activeDOFs,1)=A0;

j=1;
nADOFs=length(activeDOFs);
alphaM=0;
alphaF=-alpha;
gamma=1/2+alphaF-alphaM;
beta=0.25*((1+alphaF-alphaM)^2);

inKeq=(((1-alphaM)/(beta*(dt^2)))*M+...
    ((gamma*(1-alphaF))/(beta*dt))*C+...
    (1-alphaF)*K)\eye(nADOFs);

for t=dt:dt:time
    P1=P(t);
    P0=P(t-dt);
    Peq=P1(activeDOFs)*(1-alphaF)+P0(activeDOFs)*alphaF+...
        M*( ((1-alphaM)/(beta*dt^2))*D0+...
        ((1-alphaM)/(beta*dt))*V0-...
        ((-1+2*beta+alphaM)/(2*beta))*A0 )+...
        C*( ((gamma-gamma*alphaF)/(beta*dt))*D0-...
        ((beta-gamma+gamma*alphaF)/beta)*V0+...
        ((dt*(2*beta-gamma)*(-1+alphaF))/(2*beta))*A0 )+...
        K*( -alphaF*D0 );
    
    D1=inKeq*Peq;
    V1=(gamma/(beta*dt))*(D1-D0)+...
        (1-gamma/beta)*V0+...
        (dt-(dt*gamma)/(2*beta))*A0;
    A1=(1/(beta*dt^2))*(D1-D0)-...
        ((1/(beta*dt)))*V0+...
        (1-1/(2*beta))*A0;

    j=j+1;
    sprintf('HHT-Alpha Linear Progress: %0.2f %',(j/totalSteps)*100)
    DM(activeDOFs,j)=D1;
    VM(activeDOFs,j)=V1;
    AM(activeDOFs,j)=A1;
    D0=D1;V0=V1;A0=A1;
    
end

end