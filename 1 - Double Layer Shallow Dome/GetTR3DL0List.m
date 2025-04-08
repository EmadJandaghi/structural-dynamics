function L0List=GetTR3DL0List(refCoords,elConn)
% L0List.m
% Calculates Element initial length for 3D Truss Element
% By: Amir Hossein Namadchi     3/5/2016


nel=length(elConn); %number of Elements
L0List=zeros(nel,1);

for e=1:nel
    ni=elConn(e,1);  nj=elConn(e,2);  % element i,j (start & end ID)
    X=refCoords([ni,nj],1); % element initial X coordinates
    Y=refCoords([ni,nj],2); % element initial Y coordinates
    Z=refCoords([ni,nj],3); % element initial Z coordinates
    L0List(e)=sqrt(diff(X)^2+diff(Y)^2+diff(Z)^2);
end

end