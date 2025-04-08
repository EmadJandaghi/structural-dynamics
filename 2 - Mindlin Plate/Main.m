clc
clear
close
%% My Color Codes
cl1_POMEGRANATE=[242,38,19]/255;
cl2_MALACHITE=[0,230,64]/255;
cl3_DODGERBLUE=[25,181,254]/255;
cl4_CALIFORNIA=[248,148,6]/255;
cl5_RIPELEMON=[247,202,24]/255;
cl6_SEANCE=[154,18,179]/255;
cl7_RAZZMATAZZ=[219,10,91]/255;
cl8_REBECCAPURPLE=[102,51,153]/255;
cl9_SALEM=[30,130,76]/255;
cl10_JACKSONSPURPLE=[31,58,147]/255;
cl11_SHERPABLUE=[1,50,67]/255;
cl12_MONZA=[207,0,15]/255;
%% load Data
problemData;
Pz=-300;
rho=2.588e-4;
% load function -------------------------------
global P
P=MindlinPlateP(DOFConn,elConn,Pz,refCoords);

P=@(t) P;
% load function -------------------------------

% ---------------------------------------------

M=MindlinPlateMass(elConn,DOFConn,refCoords,rho,th,activeDOFs,1);
[K0,~]=MindlinPlateKtFint(0,elConn,DOFConn,...
    refCoords,Dh,Nx,Ny,Nxy,activeDOFs,1);
C=0*M;
%% Analyze Setting
dt=2.23e-5;
time=0.0015;

D0=zeros(nADOfs,1);
V0=zeros(nADOfs,1);
%% Processing

[DMLin,VMLin,AMLin]=BatheL(0.5,M,K0,C,activeDOFs,time,dt,D0,V0);
[DMPR,VMPR,AMPR,aTimePR]=MDEDNL(0.5,M,C,activeDOFs,time,dt,D0,V0);

[DMBNL,VMBNL,AMBNL,aTimeBNL]=BatheNLNR(0.5,M,C,activeDOFs,time,dt,D0,V0,...
    1000,1e-10);

%% Post-Processing
centralNode=(2*nex+1)*(ney/2)+nex+1;
DDOF=5*centralNode-2;
tt=0:dt:time;
hold on
hLinear=plot(tt*1000,-2.54*DMLin(DDOF,:),'Color','k','LineWidth',1.5);

hBathe=customPlot(tt*1000,-2.54*DMBNL(DDOF,:),3,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt*1000,-2.54*DMPR(DDOF,:),4,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);

ylabel('Displacement $$(cm)$$','interpreter', 'latex')               
xlabel('Time $$\times{10^{-3}} (s)$$','interpreter', 'latex')                                
% xlabel('Time ratio ({\itt\})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 22 6.5])
axis([0 time*1000 -0.05 0.59])
grid on              
grid minor

legend([hLinear,hBathe,hMDED],...
     'Linear',...
     'Bathe (Nonlinear)',...
     'MDED (Nonlinear)',...
     'Location','northoutside','Orientation','horizontal')               
% plot(tt/dt,-DDN0(DDOF,:),'Color',[0.790588, 0.201176, 0.])

%{
%% Simple Post
clk = clock;
nowDate=[num2str(clk(3)) '-' num2str(clk(2)) '-' num2str(clk(1))];
nowTime=[num2str(clk(4)) ':' num2str(clk(5)) ':' num2str(floor(clk(6)))];
nowTime1=[num2str(clk(4)) ',' num2str(clk(5)) ',' num2str(floor(clk(6)))];



fileID = fopen('result.txt','a');
fwrite(fileID,['*** ---------- ' nowDate ' ---- ' nowTime '---------- ***'])
fprintf(fileID,'\r\n');
fprintf(fileID,'dt=%8.6f sec in %5.3f sec\r\n',dt,time);
fprintf(fileID,'%dx%d with %d active DOFs\r\n',nex,ney,nADOfs);
fprintf(fileID,'  Bathe:    %8.2f sec\r\n',aTimeBNL);
fprintf(fileID,'  MEDED:    %8.2f sec\r\n',aTimePR);
fprintf(fileID,'  N-T-R:    %8.2f sec\r\n',min(aTimeBNL,aTimePR)/max(aTimeBNL,aTimePR));
fwrite(fileID,'___________________________________________________')
fprintf(fileID,'\r\n');
fclose(fileID);
%% Variable Saving
fileName=[nowDate,'(',nowTime1,')',num2str(nex),'x',num2str(ney)];
save(fileName)
%}