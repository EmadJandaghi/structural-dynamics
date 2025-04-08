% calling this function will load problem geometry into the main program
% refCoords : Problem Node Coordinates in reference Configuration
% E : Young's Modulus
% nn : number of nodes
% nDOFs : total number of degrees of freedom
% elConn : element ID connectivity matrix
% DOFConn : element DOF connectivity matrix
% fixedNodes : essential node numbers
% fixedDOFs : fixed degrees of freedoms
% activeDOFs : active degrees of freedoms
% nADOFs : number of active degrees of freedoms
% L0List : Initial length of elements
% nel : Number Of Elements
% s0List : initial stress in members
nex=4;ney=4;
Lx=10; Ly=10;
Nx=0;   Ny=0;   Nxy=0;
refCoords=[0.,0.,0;1.25,0.,0;2.5,0.,0;3.75,0.,0;5.,0.,0;6.25,0.,0;7.5,0.,0;8.75,0.,0;10.,0.,0;0.,2.5,0;1.25,2.5,0;2.5,2.5,0;3.75,2.5,0;5.,2.5,0;6.25,2.5,0;7.5,2.5,0;8.75,2.5,0;10.,2.5,0;0.,5.,0;1.25,5.,0;2.5,5.,0;3.75,5.,0;5.,5.,0;6.25,5.,0;7.5,5.,0;8.75,5.,0;10.,5.,0;0.,7.5,0;1.25,7.5,0;2.5,7.5,0;3.75,7.5,0;5.,7.5,0;6.25,7.5,0;7.5,7.5,0;8.75,7.5,0;10.,7.5,0;0.,10.,0;1.25,10.,0;2.5,10.,0;3.75,10.,0;5.,10.,0;6.25,10.,0;7.5,10.,0;8.75,10.,0;10.,10.,0;0.,1.25,0;2.5,1.25,0;5.,1.25,0;7.5,1.25,0;10.,1.25,0;0.,3.75,0;2.5,3.75,0;5.,3.75,0;7.5,3.75,0;10.,3.75,0;0.,6.25,0;2.5,6.25,0;5.,6.25,0;7.5,6.25,0;10.,6.25,0;0.,8.75,0;2.5,8.75,0;5.,8.75,0;7.5,8.75,0;10.,8.75,0];
E=1e7;   % lb/in2
nn=length(refCoords);
nDOFs=5*nn;
elConn=[1,2,3,47,12,11,10,46;10,11,12,52,21,20,19,51;19,20,21,57,30,29,28,56;28,29,30,62,39,38,37,61;3,4,5,48,14,13,12,47;12,13,14,53,23,22,21,52;21,22,23,58,32,31,30,57;30,31,32,63,41,40,39,62;5,6,7,49,16,15,14,48;14,15,16,54,25,24,23,53;23,24,25,59,34,33,32,58;32,33,34,64,43,42,41,63;7,8,9,50,18,17,16,49;16,17,18,55,27,26,25,54;25,26,27,60,36,35,34,59;34,35,36,65,45,44,43,64];
th=0.5;    %(Thickness) inch
nu=0.3;   % poissons ratio
alfa=(6/5);   % Intentionally 'alfa', not alpha (shear factor)
nel=length(elConn);

[fixedDOFs,activeDOFs]=MindlinPlateBC(nex,ney);


% Automated procedure:
DOFConn=[5*elConn(:,1)-4,5*elConn(:,1)-3,5*elConn(:,1)-2,5*elConn(:,1)-1,5*elConn(:,1),...
         5*elConn(:,2)-4,5*elConn(:,2)-3,5*elConn(:,2)-2,5*elConn(:,2)-1,5*elConn(:,2),...
         5*elConn(:,3)-4,5*elConn(:,3)-3,5*elConn(:,3)-2,5*elConn(:,3)-1,5*elConn(:,3),...
         5*elConn(:,4)-4,5*elConn(:,4)-3,5*elConn(:,4)-2,5*elConn(:,4)-1,5*elConn(:,4),...
         5*elConn(:,5)-4,5*elConn(:,5)-3,5*elConn(:,5)-2,5*elConn(:,5)-1,5*elConn(:,5),...
         5*elConn(:,6)-4,5*elConn(:,6)-3,5*elConn(:,6)-2,5*elConn(:,6)-1,5*elConn(:,6),...
         5*elConn(:,7)-4,5*elConn(:,7)-3,5*elConn(:,7)-2,5*elConn(:,7)-1,5*elConn(:,7),...
         5*elConn(:,8)-4,5*elConn(:,8)-3,5*elConn(:,8)-2,5*elConn(:,8)-1,5*elConn(:,8)];

nADOfs=length(activeDOFs);
nn=length(refCoords);

Dbar=[1,nu,0;nu,1,0;0,0,0.5*(1-nu)];
G=0.5*(E/(1+nu));
Dhp=(E*th/(1-nu^2))*Dbar;
Dhb=(E*(th^3))/(12*(1-nu^2))*Dbar;
Dhs=(G*th/alfa)*eye(2);
Dh=[Dhp,zeros(3),zeros(3,2);zeros(3),Dhb,zeros(3,2);zeros(2,6),Dhs];