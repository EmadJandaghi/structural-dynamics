clc
clear
close all
%% Memory Problem
yOn=input('Fuckin memory leak is gonna destroy your ass, continue (y/n)? :','s');
if strcmp(yOn,'n')
   return 
end
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
% (2012) - Klaus-Jürgen Bathe, Gunwoo Noh - Insight into an implicit time
% integration scheme for structural dynamics
%% Exact Solution (Modal Superposition + Static Correction)
% These soultions have been obtained based on the values reported in the
% original paper. 
D2Ex=@(t) (2.72727*(1e-7))*sin(1.0*t)+sin(1.2*t);
D3Ex=@(t) (2.72727)*sin(1.0*t)-(2.27273)*sin(1.2*t);
V2Ex=@(t) (2.72727*(1e-7))*cos(1.0*t)+1.2*cos(1.2*t);
V3Ex=@(t) (2.72727)*cos(1.0*t)-2.72727*cos(1.2*t);
A2Ex=@(t) (2.72727*(1e-7))*sin(1.0*t)-1.44*sin(1.2*t);
A3Ex=@(t) (-2.72727)*sin(1.0*t)+(3.27273*sin(1.2*t));
%% Problem Definition
global k1 wp
k1=1e7; k2=1;
m1=0;   m2=1;   m3=1;
wp=1.2;
time=10;
dt=0.2618;
activeDOFs=[2,3];
%% Initialization
K=[k1,-k1,0;
   -k1,k1+k2,-k2;
   0,-k2,k2];
M=diag([m1,m2,m3]);
C=0.0*M;
nDOFs=length(M);
[D0,V0]=deal(zeros(nDOFs,1));

M=M(activeDOFs,activeDOFs);
K=K(activeDOFs,activeDOFs);
C=C(activeDOFs,activeDOFs);
D0=D0(activeDOFs);
V0=V0(activeDOFs);
%% Processing
[DMBathe,VMBathe,AMBathe]=BatheL(0.5,M,K,C,activeDOFs,time,dt,D0,V0);
[DMMDED,VMMDED,AMMDED]=MDEDL(0.5,M,K,C,activeDOFs,time,dt,D0,V0);
[DMSSEA0,VMSSEA0,AMSSEA0]=SSEAlphaL(0.25,M,K,C,activeDOFs,time,dt,D0,V0);
[DMSSEA1,VMSSEA1,AMSSEA1]=SSEAlphaL(0.75,M,K,C,activeDOFs,time,dt,D0,V0);
[DMSSDM,VMSSDM,AMSSDM]=SSDML(-1/3,1,M,K,C,activeDOFs,time,dt,D0,V0);
%% Post-Processing
tt=0:dt:time;
ttEx=0:0.01:time;

%{
%% 1-Displacements ---------------------------------------------------------
% D2 ------------
fig1=figure;
hold on
hExact=plot(ttEx,D2Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,DMBathe(2,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,DMMDED(2,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,DMSSEA0(2,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,DMSSEA1(2,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,DMSSDM(2,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$d_{2}$$','interpreter', 'latex')               
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time min(1.2*D2Ex(ttEx)) 1.6*max(D2Ex(ttEx))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Orientation','vertical','Position',[395 168 0.1 0.2])
%% D3 ------------
fig2=figure;
hold on
hExact=plot(ttEx,D3Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,DMBathe(3,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,DMMDED(3,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,DMSSEA0(3,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,DMSSEA1(3,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,DMSSDM(3,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$d_{3}$$','interpreter', 'latex')               
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time min(1.2*D3Ex(ttEx)) 1.2*max(D3Ex(ttEx))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Orientation','vertical','Position',[252 168 0.1 0.2])
%}



%% 2-Velocities ---------------------------------------------------------
% V2 ------------
fig3=figure;
hold on
hExact=plot(ttEx,V2Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,VMBathe(2,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,VMMDED(2,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,VMSSEA0(2,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,VMSSEA1(2,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,VMSSDM(2,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$\dot{d_{2}}$$','interpreter', 'latex')                 
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time min(1.2*V2Ex(ttEx)) 1.6*max(V2Ex(ttEx))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[395 168 0.1 0.2])
%% V3 ------------
fig4=figure;
hold on
hExact=plot(ttEx,V3Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,VMBathe(3,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,VMMDED(3,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,VMSSEA0(3,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,VMSSEA1(3,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,VMSSDM(3,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$\dot{d_{3}}$$','interpreter', 'latex')            
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time min(1.2*V3Ex(ttEx)) 1.2*max(V3Ex(ttEx))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[310 77 0.1 0.2])
 
%% V2 (Different View) ------------
fig5=figure;
hold on
hExact=plot(ttEx,(1e-6)*V2Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,(1e-6)*VMBathe(2,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,(1e-6)*VMMDED(2,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,(1e-6)*VMSSEA0(2,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,(1e-6)*VMSSEA1(2,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,(1e-6)*VMSSDM(2,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$\dot{d_{2}} \times{10^{6}} $$','interpreter', 'latex')                 
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time -0.6 0.87])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[395 168 0.1 0.2])
%}


%% 3-Accelerations ---------------------------------------------------------
% A2 ------------
fig6=figure;
hold on
hExact=plot(ttEx,A2Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,AMBathe(2,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,AMMDED(2,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,AMSSEA0(2,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,AMSSEA1(2,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,AMSSDM(2,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 
ylabel('$$\ddot{d_{2}}$$','interpreter', 'latex')                
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time -24 25])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[395 168 0.1 0.2])
%% A3 ------------
fig7=figure;
hold on
hExact=plot(ttEx,A3Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,AMBathe(3,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,AMMDED(3,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,AMSSEA0(3,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,AMSSEA1(3,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,AMSSDM(3,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$\ddot{d_{3}}$$','interpreter', 'latex')            
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])

axis([0 time min(1.2*A3Ex(ttEx)) 1.2*max(A3Ex(ttEx))])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[395 168 0.1 0.2])
 
  
%% A2 (Different View) ------------
fig8=figure;
hold on
hExact=plot(ttEx,(1e-6)*A2Ex(ttEx),'Color','k','LineWidth',1.5);
 
hBathe=customPlot(tt,(1e-6)*AMBathe(2,:),1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,(1e-6)*AMMDED(2,:),2,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);
hSSEA0=customPlot(tt,(1e-6)*AMSSEA0(2,:),2,0.5*cl9_SALEM,...
                   '-.',0.5*cl9_SALEM,'d',1.1);
hSSEA1=customPlot(tt,(1e-6)*AMSSEA1(2,:),3,cl9_SALEM,...
                   '-',cl9_SALEM,'v',0.5);                
hSSDM=customPlot(tt,(1e-6)*AMSSDM(2,:),2,cl4_CALIFORNIA,...
                   '-',cl4_CALIFORNIA,'^',0.5); 

ylabel('$$\ddot{d_{2}} \times{10^{6}} $$','interpreter', 'latex')                 
xlabel('Time ({\its})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time -4.5 6])
grid on              
grid minor

legend([hExact,hBathe,hMDED,hSSEA0,hSSEA1,hSSDM],...
     'Reference',...
     'Bathe',...
     'MDED',...
     'SSE-\alpha(\mu=0.25)','SSE-\alpha(\mu=0.75)',...
     'SSDM(\alpha=-1/3)','Location','northoutside','Orientation','vertical','Position',[395 168 0.1 0.2])
%}