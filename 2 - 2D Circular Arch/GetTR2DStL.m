function K=GetTR2DStL(EList,refCoords,elConn,DOFConn,L0List,...
    areaList,activeDOFs)
%% Info
% This function constructs global stiffness matrix for a 2D Truss.
% (ref: "Structural Dynamics: Theory and Computation Book by Mario Paz")
% The code is mainly developed with the aim of increased readability, not
% efficiency. Sparse matrices might be used when memory limitations
% encountered.

% EList: Modulus of elasticity for each element (nel x 2)
% refCoords: Reference coordinate matrix of each node (nn x 2)
% elConn: Element connectivity matrix (nel x 2)
% DOFConn: DOF connectivity matrix (nel x 4)
% L0List: Length of each element (nel x 1)
% areaList: Cross-section area of each element (nel x 1)
% activeDOFs: Active degrees of freedom (nADOFs x 1) 



% Written by: A. H. Namadchi (9/2017)
%% Core
nel=size(DOFConn,1);
nn=size(refCoords,1);
nDOFs=nn*2;
K=zeros(nDOFs,nDOFs);

for e=1:nel
    DOFID=DOFConn(e,:);
    n1=elConn(e,1);n2=elConn(e,2);
    x1=refCoords(n1,1);y1=refCoords(n1,2);
    x2=refCoords(n2,1);y2=refCoords(n2,2);
    L0=L0List(e);
    c=(x2-x1)/L0;   s=(y2-y1)/L0;
    A=areaList(e);
    E=EList(e);
    Ke=(A*E/L0)*[c^2,c*s,-c^2,-c*s;
                 c*s,s^2,-c*s,-s^2;
                 -c^2,-c*s,c^2,c*s;
                 -c*s,-s^2,c*s,s^2];
    K(DOFID,DOFID)=K(DOFID,DOFID)+Ke;
    
end
K=K(activeDOFs,activeDOFs);



end