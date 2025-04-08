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
load 10-5-2018(2,7,19)10x10
%% Post-Processing
centralNode=(2*nex+1)*(ney/2)+nex+1;
DDOF=5*centralNode-2;
tt=0:dt:time;
hold on

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
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 time*1000 -0.05 0.59])
grid on              
grid minor

legend([hBathe,hMDED],...
     'Bathe',...
     'MDED',...
     'Location','northwest','Orientation','vertical','Position',[395 168 0.1 0.2])               
% plot(tt/dt,-DDN0(DDOF,:),'Color',[0.790588, 0.201176, 0.])
