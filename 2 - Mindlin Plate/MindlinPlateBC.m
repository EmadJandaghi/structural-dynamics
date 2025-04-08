function [fixedDOFs,activeDOFs]=MindlinPlateBC(nex,ney)

% u ---> 5*j-4
% v ---> 5*j-3
% w ---> 5*j-2
% thx ---> 5*j-1
% thy ---> 5*j
nnfx=2*nex+1;
nnfy=2*ney+1;
nndx=nex+1;
nndy=ney+1;
nn=nnfx*nndy+nndx*ney;
edgeBot=1:nnfx;
edgeTop=flip(nnfx*nndy:-1:nnfx*nndy-nnfx+1);
edgeLeft=[1:nnfx:nnfx*nndy-nnfx+1,nnfx*nndy+1:nndx:nn];
edgeRight=[nnfx:nnfx:nnfx*nndy,nnfx*nndy+nndx:nndx:nn];

%% All Simply Supporter
BCBot=[5*edgeBot-4,5*edgeBot-3,5*edgeBot-2,5*edgeBot-1];
BCTop=[5*edgeTop-4,5*edgeTop-3,5*edgeTop-2,5*edgeTop-1];
BCLeft=[5*edgeLeft-4,5*edgeLeft-3,5*edgeLeft-2,5*edgeLeft];
BCRight=[5*edgeRight-4,5*edgeRight-3,5*edgeRight-2,5*edgeRight];
fixedDOFs=union([BCBot,BCTop],[BCLeft,BCRight]);
activeDOFs=setdiff(1:5*nn,fixedDOFs);

end