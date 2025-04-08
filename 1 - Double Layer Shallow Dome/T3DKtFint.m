function [Kt,fint]=T3DKtFint(nD,elConn,DOFConn,refCoords,areaList...
    ,L0List,s0List,activeDOFs,E,applyBC)
% T3DKtFint.m
% Calculates Global Tangent Stiffness Matrix and
% Global internal force vector for 3D Truss Elemnt
% By: Amir Hossein Namadchi     3/5/2016


nel=length(elConn); %number of Elements
nDOFs=3*length(refCoords);  %number of Degrees of Freedom
D1=zeros(nDOFs,1);  % initial value of displacement vector
fint=zeros(nDOFs,1);    % initial value of internal force vector
Kt=sparse(nDOFs,nDOFs);    % initial value of Stiffness Matrix
%% Check nargin
if nargin<10
    applyBC=1;
end
%% nodal Displacement size Check (active or all)
if length(nD)~=nDOFs
    D1(activeDOFs)=nD;
else
    D1=nD;
end
%%  Finite Element Loop
for e=1:nel
    DOFID=DOFConn(e,:); % element DOF identifiers
    LL0=L0List(e)^2;  % initial length squared
    A=areaList(e);  %   element cross-section area
    s0=s0List(e);    %   element initial Stress
    n1=elConn(e,1);  n2=elConn(e,2);  % element i,j (start & end ID)
    
    %--------% element initial coordinates----
    X1=refCoords(n1,1);X2=refCoords(n2,1);
    Y1=refCoords(n1,2);Y2=refCoords(n2,2);
    Z1=refCoords(n1,3);Z2=refCoords(n2,3);
    %--------% element nodal deflection-------
    ux1=D1(3*n1-2);ux2=D1(3*n2-2);
    vy1=D1(3*n1-1);vy2=D1(3*n2-1);
    wz1=D1(3*n1);wz2=D1(3*n2);
    %--------% element current coordinates----
    x1=X1+ux1;  x2=X2+ux2;
    y1=Y1+vy1;  y2=Y2+vy2;
    z1=Z1+wz1;  z2=Z2+wz2;
    % ----------------------------------------
    LL=(x2-x1)^2+(y2-y1)^2+(z2-z1)^2;
    strain=(LL-LL0)/(2*LL0);    % Green Strain
    B=[-(x2-x1);-(y2-y1);-(z2-z1);(x2-x1);(y2-y1);(z2-z1)];
    BB=(B*B.')/LL0;
    J=[1,0,0,-1,0,0;...
        0,1,0,0,-1,0;...
        0,0,1,0,0,-1;...
        -1,0,0,1,0,0;...
        0,-1,0,0,1,0;...
        0,0,-1,0,0,1];
    f0=s0*A;
    f=f0+E*A*strain;
    finte=(f*B)/L0List(e);  % element internal force vector
    Kte=(A*E/L0List(e))*BB+...
        (f/L0List(e))*J;  % element tangent stiffness Matrix
    Kt(DOFID,DOFID)=Kt(DOFID,DOFID)+Kte;    % Global tangent stiffness Matrix
    fint(DOFID)=fint(DOFID)+finte;  % global internal force vector
end
if applyBC==1
    Kt=Kt(activeDOFs,activeDOFs);
    fint=fint(activeDOFs);
end
end