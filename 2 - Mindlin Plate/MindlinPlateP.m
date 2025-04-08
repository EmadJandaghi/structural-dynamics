function P=MindlinPlateP(DOFConn,elConn,Pz,refCoords)
nel=length(elConn);
nn=length(refCoords);
nDOFs=5*nn;
P=zeros(nDOFs,1);

for e=1:nel
 elID=elConn(e,:);
 DOFID=DOFConn(e,:);
 elCoords=refCoords(elID,[1,2]);
 pe=GetElementP(Pz,elCoords);
 P(DOFID)=P(DOFID)+pe;
end
end

function pe=GetElementP(Pz,elCoords)
wg=[2,0,0;1,1,0;5/9,8/9,5/9];
xg=[0,0,0;-(1/sqrt(3)),1/sqrt(3),0;-sqrt(3/5),0,sqrt(3/5)];
ngp=3;
% Initializing
pe= zeros(40,1);

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
      dNVxi=[-((-1+eta)*(eta+2*xi))/4,(-1+eta)*xi,...
          ((-1+eta)*(eta-2*xi))/4,(1-eta^2)/2,((1+eta)*(eta+2*xi))/4,...
          (-1-eta)*xi,-((1+eta)*(eta-2*xi))/4,(-1+eta^2)/2];
      dNVeta=[-((-1+xi)*(2*eta+xi))/4,(-1+xi^2)/2,...
          ((2*eta-xi)*(1+xi))/4,-(eta*(1+xi)),((1+xi)*(2*eta+xi))/4,...
          (1-xi^2)/2,-((2*eta-xi)*(-1+xi))/4,eta*(-1+xi)];
      J=[dNVxi;dNVeta]*elCoords;
      NM=[NV(1),0,0,0,0,NV(2),0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),...
          0,0,0,0,NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8),0,0,0,0;0,NV(1),0,0,...
          0,0,NV(2),0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,...
          NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8),0,0,0;0,0,NV(1),0,0,0,0,NV(2),...
          0,0,0,0,NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,...
          0,NV(7),0,0,0,0,NV(8),0,0;0,0,0,NV(1),0,0,0,0,NV(2),0,0,0,0,...
          NV(3),0,0,0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,0,NV(7),...
          0,0,0,0,NV(8),0;0,0,0,0,NV(1),0,0,0,0,NV(2),0,0,0,0,NV(3),0,0,...
          0,0,NV(4),0,0,0,0,NV(5),0,0,0,0,NV(6),0,0,0,0,NV(7),0,0,0,0,NV(8)];
      pe=pe+wi*wj*det(J)*(transpose(NM)*[0;0;Pz;0;0]);
  end
end

end