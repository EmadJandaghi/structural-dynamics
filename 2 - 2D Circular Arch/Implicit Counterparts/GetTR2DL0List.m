function L0List=GetTR2DL0List(refCoords,elConn)
%% Info
% This function calculates the length of each elements for a 2D Truss.
% The code is mainly developed with the aim of increased readability, not
% efficiency.

% refCoords: Reference coordinate matrix of each node (nn x 2)
% elConn: Element connectivity matrix (nel x 2)


% Written by: A. H. Namadchi (9/2017)
%% Core
nel=size(elConn,1);
L0List=zeros(nel,1);
for e=1:nel
   i=elConn(e,1); j=elConn(e,2);
   X1=refCoords(i,1);X2=refCoords(j,1);
   Y1=refCoords(i,2);Y2=refCoords(j,2);
   L0List(e)=sqrt((X2-X1)^2+(Y2-Y1)^2);
   
end


end