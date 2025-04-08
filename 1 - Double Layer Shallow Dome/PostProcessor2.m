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
load 10-5-2018(1,4,16)
% load function -------------------------------

% ---------------------------------------------

M=GetTR3DMass(rho,refCoords,DOFConn,L0List,areaList,activeDOFs);
C=0*M;
%% Analyze Setting
dt=0.00025;
time=0.2;

D0=zeros(nADOfs,1);
V0=zeros(nADOfs,1);


%% Post-Processing

DDOF=3*362;
tt=0:dt:time;
hold on


hBathe=customPlot(tt,100*DMBNL(DDOF,:),23,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(tt,100*DMPR(DDOF,:),25,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);

ylabel('Displacement $$(cm)$$','interpreter', 'latex')               
xlabel('Time $$(s)$$','interpreter', 'latex')                                
% xlabel('Time ratio ({\itt\})')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12.5)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time -35 35])
grid on              
grid minor

legend([hBathe,hMDED],...     
     ['Bathe   - CPU Time: ',num2str(ceil(aTimeBNL)),' s'],...
     ['MDED - CPU Time: ',num2str(ceil(aTimePR)),' s'],...
     'Location','northeast','Orientation','vertical','Position',[328 195 0.1 0.2])               
% plot(tt/dt,-DDN0(DDOF,:),'Color',[0.790588, 0.201176, 0.])
