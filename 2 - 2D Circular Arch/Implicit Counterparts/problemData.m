%% Calling this function will load problem geometry into the main program

% refCoords: Reference coordinate matrix of each node (nn x 2)
% EList : Young's Modulus List (nel x 1)
% nn : number of nodes
% nDOFs : total number of degrees of freedom
% elConn: Element connectivity matrix (nel x 2)
% DOFConn: DOF connectivity matrix (nel x 4)
% fixedNodes : fixed nodes number list
% fixedDOFs : fixed degrees of freedoms list
% activeDOFs : active degrees of freedoms (nADOFs x 1)
% nADOFs : number of active degrees of freedoms
% L0List : Initial length of elements
% nel : Number Of Elements

global nDOFs

refCoords=[-33.9411,33.9411;-31.1735,36.4995;-28.2137,38.8328;
    -25.0799,40.9267;-21.7915,42.7683;-18.3688,44.3462;-14.8328,45.6507;
    -11.2054,46.6738;-7.50885,47.409;-3.76604,47.852;0.,48.;3.76604,47.852;
    7.50885,47.409;11.2054,46.6738;14.8328,45.6507;18.3688,44.3462;
    21.7915,42.7683;25.0799,40.9267;28.2137,38.8328;31.1735,36.4995;
    33.9411,33.9411;-35.3553,35.3553;-32.4724,38.0203;-29.3893,40.4508;
    -26.1249,42.632;-22.6995,44.5503;-19.1342,46.194;-15.4508,47.5528;
    -11.6723,48.6185;-7.82172,49.3844;-3.92295,49.8459;0.,50.;
    3.92295,49.8459;7.82172,49.3844;11.6723,48.6185;15.4508,47.5528;
    19.1342,46.194;22.6995,44.5503;26.1249,42.632;29.3893,40.4508;
    32.4724,38.0203;35.3553,35.3553]*(0.01);   % m
       
elConn=[1,22;2,23;3,24;4,25;5,26;6,27;7,28;8,29;9,30;10,31;11,32;12,33;
    13,34;14,35;15,36;16,37;17,38;18,39;19,40;20,41;21,42;1,2;2,3;3,4;
    4,5;5,6;6,7;7,8;8,9;9,10;10,11;11,12;12,13;13,14;14,15;15,16;16,17;
    17,18;18,19;19,20;20,21;22,23;23,24;24,25;25,26;26,27;27,28;28,29;
    29,30;30,31;31,32;32,33;33,34;34,35;35,36;36,37;37,38;38,39;39,40;
    40,41;41,42;1,23;2,24;3,25;4,26;5,27;6,28;7,29;8,30;9,31;10,32;
    11,33;12,34;13,35;14,36;15,37;16,38;17,39;18,40;19,41;20,42;22,2;23,3;
    24,4;25,5;26,6;27,7;28,8;29,9;30,10;31,11;32,12;33,13;34,14;35,15;
    36,16;37,17;38,18;39,19;40,20;41,21];

E=200e9;   % N/m2
rho=8000;    % kg/m3
A=0.00025;    % m2


nn=length(refCoords);
nDOFs=2*nn;
nel=length(elConn);

areaList=A*ones(nel,1); % in2
EList=E*ones(nel,1);
fixedNodes=[1,21];
fixedDOFs=[2*fixedNodes-1,2*fixedNodes];


% Automated procedure:
DOFConn=[2*elConn(:,1)-1,2*elConn(:,1),...
         2*elConn(:,2)-1,2*elConn(:,2)];
L0List=GetTR2DL0List(refCoords,elConn);

activeDOFs=setdiff(1:nDOFs,fixedDOFs);
nADOfs=length(activeDOFs);


