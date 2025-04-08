function [DM,VM,AM]=SSEAlphaL(mu,M,K,C,activeDOFs,time,dt,D0,V0)
%% Info
% Based on:
% (2015) - C. Kolay, James M. Ricles - Assessment of explicit and
% semi-explicit classes of model-based algorithms for direct integration in
% structural dynamics


% In this method, the displacement difference equation is obtained by
% modifying that of the NE method with two model-based integration
% parameter matrices, whereas velocity equation of the N-Beta method is
% adopted for velocity. The controllable numerical dissipation is
% introduced through the weighted equations of motion of the KR-alpha
% method. Hereafter, this subfamily of the SE-alpha method is
% referred to as the single-parameter semi-explicit-alpha (SSE-alpha)
% method.

% mu=0: Maximum High frequency dissipation
% mu=1: No Dissipation



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
% Preliminaries
j=1;
nADOFs=length(activeDOFs);
IM=eye(nADOFs);


% Weighting Factors
alphaM=(2*mu-1)/(mu+1);
alphaF=mu/(mu+1);
gamma=0.5-alphaM+alphaF;
beta=0.25*(0.5+gamma)^2;

% Integration Matrices
alpha=M+gamma*dt*C+beta*(dt^2)*K;
delta=(gamma-2*beta);
alphaHat=alphaM*M+alphaM*gamma*dt*C+alphaF*beta*(dt^2)*K;
alpha1=alpha\(M+gamma*dt*C);
alpha2=0.5*(alpha\(M+delta*dt*C));
alpha3=alpha\alphaHat;

% Factorizations
% I didn't use inv(A) since it's less accurate
% inCeq=((1/(dt*gamma))*(M*(IM-alpha3))+(1-alphaF)*C)\IM;
inMeq=(M*(IM-alpha3)+dt*gamma*(1-alphaF)*C)\IM;


for t=dt:dt:time
    P1=P(t);
    P0=P(t-dt);
    
    D1=D0+dt*(alpha1*V0)+(dt^2)*(alpha2*A0);
    Df=(1-alphaF)*D1+alphaF*D0;
    Pf=(1-alphaF)*P1(activeDOFs)+alphaF*P0(activeDOFs);
    
    A1=inMeq*(Pf-K*Df-M*alpha3*A0-C*(V0+dt*((1-alphaF)*(1-gamma)*A0)));
    V1=V0+dt*((1-gamma)*A0+gamma*A1);
    
    j=j+1;
    sprintf('SSE-Alpha Linear Progress: %0.2f %',(j/totalSteps)*100)
    DM(activeDOFs,j)=D1;
    VM(activeDOFs,j)=V1;
    AM(activeDOFs,j)=A1;
    D0=D1;V0=V1;A0=A1;
    
end

end