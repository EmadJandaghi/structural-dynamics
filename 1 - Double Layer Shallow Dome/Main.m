clc
clear
%% load Data
problemData;

% load function -------------------------------

% ---------------------------------------------
[Kt,~]=T3DKtFint(0,elConn,DOFConn,refCoords,areaList...
    ,L0List,s0List,activeDOFs,E,1);
M=GetTR3DMass(rho,refCoords,DOFConn,L0List,areaList,activeDOFs);
C=0*M;
%% Analyze Setting
dt=0.00025;
time=0.2;

D0=zeros(nADOfs,1);
V0=zeros(nADOfs,1);
tolR=1e-6;
tolK=1e-12;
maxDRIteration=6000;
%% Processing


[DMPR,VMPR,AMPR,aTimePR]=MDEDNL(0.5,M,C,activeDOFs,time,dt,D0,V0);

[DMBNL,VMBNL,AMBNL,aTimeBNL]=BatheNLNR(0.5,M,C,activeDOFs,time,dt,D0,V0,...
    1000,1e-10);

%% Post-Processing
DDOF=3*362;
tt=0:dt:time;
hold on
plot(tt,-DMPR(DDOF,:),'-')
plot(tt,-DMBNL(DDOF,:),'.-')
legend('MDED','Bathe')
axis([0 time 1.1*min(-DMBNL(DDOF,:)) 1.1*max(-DMBNL(DDOF,:))])
% plot(tt/dt,-DDN0(DDOF,:),'Color',[0.790588, 0.201176, 0.])
%% Simple Post
clk = clock;
nowDate=[num2str(clk(3)) '-' num2str(clk(2)) '-' num2str(clk(1))];
nowTime=[num2str(clk(4)) ':' num2str(clk(5)) ':' num2str(floor(clk(6)))];
nowTime1=[num2str(clk(4)) ',' num2str(clk(5)) ',' num2str(floor(clk(6)))];



fileID = fopen('result.txt','a');
fwrite(fileID,['*** ---------- ' nowDate ' ---- ' nowTime '---------- ***'])
fprintf(fileID,'\r\n');
fprintf(fileID,'dt=%8.6f sec in %5.3f sec',dt,time);
fprintf(fileID,'\r\n');
fprintf(fileID,'  Bathe:    %8.2f sec\r\n',aTimeBNL);
fprintf(fileID,'  MEDED:    %8.2f sec\r\n',aTimePR);
fprintf(fileID,'  N-T-R:    %8.2f sec\r\n',min(aTimeBNL,aTimePR)/max(aTimeBNL,aTimePR));
fwrite(fileID,'___________________________________________________')
fprintf(fileID,'\r\n');
fclose(fileID);
%% Variable Saving
fileName=[nowDate,'(',nowTime1,')'];
save(fileName)
