function M=GetTR3DMass(rho,refCoords,DOFConn,L0List,areaList,activeDOFs)
%% Info
% This function constructs consistent global mass matrix for a 3D Truss.
% (ref: "Structural Dynamics: Theory and Computation Book by Mario Paz")
% The code is mainly developed with the aim of increased readability, not
% efficiency. Sparse matrices might be used when memory limitations
% encountered.

% rho: Mass density
% refCoords: Reference coordinate matrix of each node (nn x 3)
% DOFConn: DOF connectivity matrix (nel x 6)
% L0List: Length of each element (nel x 1)
% areaList: Cross-section area of each element (nel x 1)
% activeDOFs: Active degrees of freedom (nADOFs x 1) 



% Written by: A. H. Namadchi (9/2017)
%% Core
nel=size(DOFConn,1);
nn=size(refCoords,1);
nDOFs=nn*3;
M=sparse(nDOFs,nDOFs);

for e=1:nel
    DOFID=DOFConn(e,:);
    L0=L0List(e);
    A=areaList(e);
    Me=(1/6)*(rho*A*L0)*[2,0,0,1,0,0;
                        0,2,0,0,1,0;
                        0,0,2,0,0,1;
                        1,0,0,2,0,0;
                        0,1,0,0,2,0;
                        0,0,1,0,0,2];
    M(DOFID,DOFID)=M(DOFID,DOFID)+Me;
    
end
M=M(activeDOFs,activeDOFs);



end