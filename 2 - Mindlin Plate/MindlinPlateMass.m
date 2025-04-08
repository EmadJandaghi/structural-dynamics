function M=MindlinPlateMass(elConn,DOFConn,refCoords,rho,th,activeDOFs,applyBC)
% MindlinPlateMass.m
% Calculates Global Mass Matrix of Mindlin Plate
% I've found nothing about mindlin geometrically nonlinear mindlin mass
% matrix. However, I choose to follow the similar procedure which is used
% in this Book on page 300:
% -------------------------------------------------------
% 26 - Advances in Computational Dynamics of Particles, Materials and
% Structures - Jason Har, Kumar K. Tamma (2012) 
% -------------------------------------------------------
% This is a finite element analysis of the geometrically nonlinear
% behaviour of plates using a Mindlin formulation with the assumption of small rotations.
% Based On:----------------------------------------
% (1980) - A. Pica, R.D. Wood, E. Hinton - 
% Finite element analysis of geometrically nonlinear plate behaviour
% using a mindlin formulation
% -------------------------------------------------
% Here, 8 Node serendipity Element (QS) is employed

nel=length(elConn); %number of Elements
nDOFs=5*length(refCoords);  %number of Degrees of Freedom
M=sparse(nDOFs,nDOFs);    % initial value of Mass Matrix
%% Check nargin
if nargin<7
    applyBC=1;
end
%%  Finite Element Loop
for e=1:nel
   elID=elConn(e,:);
   DOFID=DOFConn(e,:);
   elCoords=refCoords(elID,[1,2]);
   me=GetElementM(elCoords,rho,th);
   M(DOFID,DOFID)=M(DOFID,DOFID)+me;
end
if applyBC==1
    M=M(activeDOFs,activeDOFs);
end
end

function me=GetElementM(elCoords,rho,th)
wg=[2,0,0;1,1,0;5/9,8/9,5/9];
xg=[0,0,0;-(1/sqrt(3)),1/sqrt(3),0;-sqrt(3/5),0,sqrt(3/5)];
ngp=3;
% Initializing
me=zeros(40);

for i=1:ngp
    wi=wg(ngp,i);
    xi=xg(ngp,i);
  for j=1:ngp
      wj=wg(ngp,j);
      eta=xg(ngp,j);
      NV=[-((-1+eta)*(-1+xi)*(1+eta+xi))/4,...
          ((-1+eta)*(-1+xi^2))/2,((-1+eta)*(1+eta-xi)*(1+xi))/4,...
          -((-1+eta^2)*(1+xi))/2,((1+eta)*(1+xi)*(-1+eta+xi))/4,...
          -((1+eta)*(-1+xi^2))/2,((1+eta)*(-1+xi)*(1-eta+xi))/4,...
          ((-1+eta^2)*(-1+xi))/2];
      NM=[NV(1),0,0,0,0,NV(2),0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),...
          0,0,0,0,NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8),0,0,0,0;0,NV(1),0,0,...
          0,0,NV(2),0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,...
          NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8),0,0,0;0,0,NV(1),0,0,0,0,NV(2),...
          0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,...
          0,NV(7),0,0,0,0,NV(8),0,0;0,0,0,NV(1),0,0,0,0,NV(2),0,0,0,0,...
          NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,0,NV(7),...
          0,0,0,0,NV(8),0;0,0,0,0,NV(1),0,0,0,0,NV(2),0,0,0,0,NV(3),0,0,...
          0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8)];
      
      dNVxi=[-((-1+eta)*(eta+2*xi))/4,(-1+eta)*xi,...
          ((-1+eta)*(eta-2*xi))/4,(1-eta^2)/2,((1+eta)*(eta+2*xi))/4,...
          (-1-eta)*xi,-((1+eta)*(eta-2*xi))/4,(-1+eta^2)/2];
      dNVeta=[-((-1+xi)*(2*eta+xi))/4,(-1+xi^2)/2,...
          ((2*eta-xi)*(1+xi))/4,-(eta*(1+xi)),((1+xi)*(2*eta+xi))/4,...
          (1-xi^2)/2,-((2*eta-xi)*(-1+xi))/4,eta*(-1+xi)];
      J=[dNVxi;dNVeta]*elCoords;
      
      
      me=me+wi*wj*det(J)*(transpose(NM)*...
          [rho*th,0,0,0,0;
           0,rho*th,0,0,0;
           0,0,rho*th,0,0;
           0,0,0,(rho*th^3)/12,0;
           0,0,0,0,(rho*th^3)/12]*NM);
  end
end

end