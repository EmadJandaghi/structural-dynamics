clc
clear
close all
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
%% Reference
% (2009) - Huu-Tai Thai, Seung-Eock Kim - Large deflection inelastic
% analysis of space trusses using generalized displacement control method

% 5.1. Two-dimensional circular arch truss

% ******************** My Note ********************:
% Only Geometry and Configuration have been adopted from the paper.
% ******************** ^^^^^^^ ********************
%% Problem Definition
problemData;
K=GetTR2DStL(EList,refCoords,elConn,DOFConn,L0List,areaList,activeDOFs);
M=GetTR2DMass(rho,refCoords,DOFConn,L0List,areaList,activeDOFs);
C=0*M;
wMin=sqrt(eigs(K,M,1,'sm'));wMax=sqrt(eigs(K,M,1,'lm'));
TMin=2*pi/wMax;TMax=2*pi/wMin;
time=0.15;
dt=0.0075;

D0=0*(activeDOFs');
V0=D0*0;         % initial Conditions
%% Exact Solution (Obtained by the Bathe Method)
dtEx=0.001;
[DMEx,VMEx,AMEx]=BatheL(0.5,M,K,C,activeDOFs,time,dtEx,D0,V0);
%% Processing
[DMBathe,VMBathe,AMBathe]=BatheL(0.5,M,K,C,activeDOFs,time,dt,D0,V0);
[DMMDED,VMMDED,AMMDED]=MDEDL(0.5,M,K,C,activeDOFs,time,dt,D0,V0);
[DMGA0,VMGA0,AMGA0]=GenAlphaL(0.25,M,K,C,activeDOFs,time,dt,D0,V0);
[DMGA1,VMGA1,AMGA1]=GenAlphaL(0.75,M,K,C,activeDOFs,time,dt,D0,V0);
[DMHHT,VMHHT,AMHHT]=HHTAlphaL(-1/3,M,K,C,activeDOFs,time,dt,D0,V0);
%% Post-Processing
tt=0:dt:time;
ttEx=0:dtEx:time;
%% 1-Displacements ---------------------------------------------------------
% fig1=figure;
% hold on
% dofToPlotD=2*32;
% hExact=plot(ttEx,DMEx(dofToPlotD,:),'Color','k','LineWidth',1.5);
% 
% hBathe=customPlot(tt,DMBathe(dofToPlotD,:),1,cl1_POMEGRANATE,...
%                    '-',cl1_POMEGRANATE,'s',0.5);
% hMDED=customPlot(tt,DMMDED(dofToPlotD,:),2,cl3_DODGERBLUE,...
%                    '--',cl3_DODGERBLUE,'o',1.0);
% hGA0=customPlot(tt,DMGA0(dofToPlotD,:),1,0.5*cl9_SALEM,...
%                    '-.',0.5*cl9_SALEM,'d',1.1);
% hGA1=customPlot(tt,DMGA1(dofToPlotD,:),1,cl9_SALEM,...
%                    '-',cl9_SALEM,'v',0.5);                
% hHHT=customPlot(tt,DMHHT(dofToPlotD,:),1,cl4_CALIFORNIA,...
%                    '-',cl4_CALIFORNIA,'^',0.5); 
% 
% ylabel('Displacement $$(m)$$','interpreter', 'latex')               
% xlabel('Time ({\its})')
% set(gca,'FontName','Times New Roman')
% set(gca,'box','on')
% set(gcf, 'Units','centimeters', 'Position',[5 5 22 6.5])
% axis([0 time 0 0.026])
% grid on              
% grid minor
% 
% legend([hExact,hBathe,hMDED,hGA0,hGA1,hHHT],...
%      'Reference',...
%      'Bathe',...
%      'MDED',...
%      'G-\alpha(\mu=0.25)','G-\alpha(\mu=0.75)',...
%      'HHT(\alpha=-1/3)','Location','northoutside','Orientation','horizontal')
%  
%% 2-Velocities ---------------------------------------------------------
% Different View
% fig2=figure;
% hold on
% dofToPlotV=2*32;
% hExact=plot(ttEx,VMEx(dofToPlotV,:),'Color','k','LineWidth',1.5);
% 
% hBathe=customPlot(tt,VMBathe(dofToPlotV,:),1,cl1_POMEGRANATE,...
%                    '-',cl1_POMEGRANATE,'s',0.5);
% hMDED=customPlot(tt,VMMDED(dofToPlotV,:),2,cl3_DODGERBLUE,...
%                    '--',cl3_DODGERBLUE,'o',1.0);
% hGA0=customPlot(tt,VMGA0(dofToPlotV,:),1,0.5*cl9_SALEM,...
%                    '-.',0.5*cl9_SALEM,'d',1.1);
% hGA1=customPlot(tt,VMGA1(dofToPlotV,:),1,cl9_SALEM,...
%                    '-',cl9_SALEM,'v',0.5);                
% hHHT=customPlot(tt,VMHHT(dofToPlotV,:),1,cl4_CALIFORNIA,...
%                    '-',cl4_CALIFORNIA,'^',0.5); 
% 
% ylabel('Velocity $$(m/s)$$','interpreter', 'latex')               
% xlabel('Time ({\its})')
% set(gca,'FontName','Times New Roman')
% set(gca,'box','on')
% set(gcf, 'Units','centimeters', 'Position',[5 5 22 6.5])
% axis([0 time -2000 7000])
% grid on              
% grid minor
% 
% legend([hExact,hBathe,hMDED,hGA0,hGA1,hHHT],...
%      'Reference',...
%      'Bathe',...
%      'MDED',...
%      'G-\alpha(\mu=0.25)','G-\alpha(\mu=0.75)',...
%      'HHT(\alpha=-1/3)','Location','northoutside','Orientation','horizontal')
%  
%% Standard View
fig3=figure;
hold on
dofToPlotV=2*32;
hExact=plot(ttEx,VMEx(dofToPlotV,:),'Color','k','LineWidth',1.5);

hBathe=customPlot(tt,VMBathe(dofToPlotV,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,VMMDED(dofToPlotV,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hGA0=customPlot(tt,VMGA0(dofToPlotV,:),1,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hGA1=customPlot(tt,VMGA1(dofToPlotV,:),1,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hHHT=customPlot(tt,VMHHT(dofToPlotV,:),1,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('Velocity $$(m/s)$$','interpreter', 'latex')               
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time 1.1*min(VMEx(dofToPlotV,:)) 1.2*max(VMEx(dofToPlotV,:))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hGA0,hGA1,hHHT],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'G-\alpha(\mu=0.25)','G-\alpha(\mu=0.75)',...
     'HHT(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[358 168 0.1 0.2])
 
%% 3-Accelerations --------------------------------------------------------
% Different View
% fig4=figure;
% hold on
% dofToPlotA=2*32;
% hExact=plot(ttEx,(1e-6)*AMEx(dofToPlotA,:),'Color','k','LineWidth',1.5);
% 
% hBathe=customPlot(tt,(1e-6)*AMBathe(dofToPlotA,:),1,cl1_POMEGRANATE,...
%                    '-',cl1_POMEGRANATE,'s',0.5);
% hMDED=customPlot(tt,(1e-6)*AMMDED(dofToPlotA,:),2,cl3_DODGERBLUE,...
%                    '--',cl3_DODGERBLUE,'o',1.0);
% hGA0=customPlot(tt,(1e-6)*AMGA0(dofToPlotA,:),1,0.5*cl9_SALEM,...
%                    '-.',0.5*cl9_SALEM,'d',1.1);
% hGA1=customPlot(tt,(1e-6)*AMGA1(dofToPlotA,:),1,cl9_SALEM,...
%                    '-',cl9_SALEM,'v',0.5);                
% hHHT=customPlot(tt,(1e-6)*AMHHT(dofToPlotA,:),1,cl4_CALIFORNIA,...
%                    '-',cl4_CALIFORNIA,'^',0.5); 
% 
% ylabel('Acceleration $$\times{10^{6}} (m/s^{2})$$','interpreter', 'latex')               
% xlabel('Time ({\its})')
% set(gca,'FontName','Times New Roman')
% set(gca,'box','on')
% set(gcf, 'Units','centimeters', 'Position',[5 5 22 6.5])
% % axis([0 time -2000 7000])
% grid on              
% grid minor
% 
% legend([hExact,hBathe,hMDED,hGA0,hGA1,hHHT],...
%      'Reference',...
%      'Bathe',...
%      'MDED',...
%      'G-\alpha(\mu=0.25)','G-\alpha(\mu=0.75)',...
%      'HHT(\alpha=-1/3)','Location','northoutside','Orientation','horizontal')
%  
% Standard View
fig5=figure;
hold on
dofToPlotA=2*32;
hExact=plot(ttEx,AMEx(dofToPlotA,:),'Color','k','LineWidth',1.5);

hBathe=customPlot(tt,AMBathe(dofToPlotA,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,AMMDED(dofToPlotA,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hGA0=customPlot(tt,AMGA0(dofToPlotA,:),1,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hGA1=customPlot(tt,AMGA1(dofToPlotA,:),1,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hHHT=customPlot(tt,AMHHT(dofToPlotA,:),1,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('Acceleration $$(m/s^{2})$$','interpreter', 'latex')               
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time 1.0*min(AMEx(dofToPlotA,:)) 1.1*max(AMEx(dofToPlotA,:))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hGA0,hGA1,hHHT],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'G-\alpha(\mu=0.25)','G-\alpha(\mu=0.75)',...
     'HHT(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[358 168 0.1 0.2])