clc
clear
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
load 10-5-2018(1,4,16)
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

[DMLin,VMLin,AMLin]=BatheL(0.5,M,Kt,C,activeDOFs,time,dt,D0,V0);

%% Post-Processing

DDOF=3*362;
tt=0:dt:time;
hold on
hLinear=plot(tt,DMLin(DDOF,:),'Color','k','LineWidth',1.5);

hBathe=customPlot(tt,DMBNL(DDOF,:),23,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,DMPR(DDOF,:),25,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);

ylabel('Displacement $$(cm)$$','interpreter', 'latex')               
xlabel('Time $$(s)$$','interpreter', 'latex')                                
% xlabel('Time ratio ({\itt\})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 22 6.5])
% axis([0 time*1000 -0.05 0.59])
grid on              
grid minor

legend([hLinear,hBathe,hMDED],...
     'Linear',...
     'Bathe (Nonlinear)',...
     'MDED (Nonlinear)',...
     'Location','northoutside','Orientation','horizontal')               
% plot(tt/dt,-DDN0(DDOF,:),'Color',[0.790588, 0.201176, 0.])
