function [DM,VM,AM]=SSDML(alpha,sigma,M,K,C,activeDOFs,time,dt,D0,V0)
%% Info
% Based on:
% (2014) - S. Y. Chang - A family of noniterative integration methods with
% desired numerical dissipation


% A new family of unconditionally stable integration methods for structural
% dynamics has been developed, which possesses the favorable numerical
% dissipation properties that can be continuously controlled. In
% particular, it can have zero damping. The most important improvement of
% this family method is that it involves no nonlinear iterations for each
% time step, and thus it is very computationally efficient when compared
% with a general second-order accurate integration method, such as the
% constant average acceleration method.

% sigma: stability amplification factor (1 or 2)

% -1/3 <= alpha <= 0
% alpha=-1/3 --> more dissipative, more period error
% alpha=-1/6 --> less dissipative, less period error
% alpha=0    --> no dissipation, least period error (CAAM)

% *** In the paper, alpha=-1/3 was chosen for the numerical examples.

% Written by: A. H. Namadchi (10/2017)
%% Core
P0=P(0);
nDOFs=length(P0);
A0=M\(P0(activeDOFs)-K*D0-C*V0);
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


% Weighting Factors
beta=0.25*((1-alpha)^2);
gamma=0.5-alpha;

% Integration Matrices
beta1=IM + (M+gamma*dt*C+(1+alpha)*beta*sigma*(dt^2)*K)\...
            (alpha*beta*sigma*((dt^2)*K));
beta2=(M+gamma*dt*C+(1+alpha)*beta*sigma*(dt^2)*K)\(M+gamma*dt*C);
beta3=(M+gamma*dt*C+(1+alpha)*beta*sigma*(dt^2)*K)\...
      (0.5*M-(beta-0.5*gamma)*dt*C); 
  

% Factorizations
% I didn't use inv(A) since it's less accurate
inCeq=(C+(1/(dt*gamma))*M)\IM;

%% Starting procedure (alpha=0)
beta_0=0.25*((1-(0))^2);
gamma_0=0.5-(0);
beta1_0=IM + (M+gamma_0*dt*C+(1+(0))*beta_0*sigma*(dt^2)*K)\...
            ((0)*beta_0*sigma*((dt^2)*K));
beta2_0=(M+gamma_0*dt*C+(1+(0))*beta_0*sigma*(dt^2)*K)\(M+gamma_0*dt*C);
beta3_0=(M+gamma_0*dt*C+(1+(0))*beta_0*sigma*(dt^2)*K)\...
      (0.5*M-(beta_0-0.5*gamma_0)*dt*C);
D1=beta1_0*D0+beta2_0*dt*V0+beta3_0*(dt^2)*A0;
Da=(1+(0))*D1-(0)*D0;
P1=P(dt);P0=P(0);
Pa=(1+(0))*P1(activeDOFs)-(0)*P0(activeDOFs);
V1=(C+(1/(dt*gamma_0))*M)\...
   (Pa-K*Da-((-1+gamma_0)/gamma_0)*(M*A0)+(1/(dt*gamma_0))*(M*V0));
A1=((gamma_0-1)/(gamma_0))*A0+(1/(dt*gamma_0))*(V1-V0);
DM(activeDOFs,2)=D1;
VM(activeDOFs,2)=V1;
AM(activeDOFs,2)=A1;
j=j+1;
Dm1=D0;D0=D1;V0=V1;A0=A1;
clear beta_0  gamma_0 beta1_0 beta2_0 beta3_0


for t=2*dt:dt:time
    P1=P(t);
    P0=P(t-dt);
    D1=(IM-beta1)*Dm1+beta1*D0+beta2*dt*V0+beta3*(dt^2)*A0;
    Da=(1+alpha)*D1-(alpha)*D0;    
    Pa=(1+alpha)*P1(activeDOFs)-(alpha)*P0(activeDOFs);
    V1=inCeq*(Pa-K*Da-((-1+gamma)/gamma)*(M*A0)+(1/(dt*gamma))*(M*V0));
    A1=((gamma-1)/(gamma))*A0+(1/(dt*gamma))*(V1-V0);    
    
    j=j+1;
    sprintf('SSDM Linear Progress: %0.2f %',(j/totalSteps)*100)
    DM(activeDOFs,j)=D1;
    VM(activeDOFs,j)=V1;
    AM(activeDOFs,j)=A1;
    Dm1=D0;D0=D1;V0=V1;A0=A1;
    
end

end