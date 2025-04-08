clc
clear
close
%% Simple


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

%% Post-Processing
nDOFs=[0,193,469,865,1381,2017,2773,3649,4645,5761];
tBathe=[0,7.94,19.43,37.96,65.42,109.50,178.63,282.69,428.69,622.27];
tMDED=[0,0.57,2.14,5.89,13.20,27.06,48.35,81.08,129.39,205.41];

fig1=figure;
hold on

hBathe=customPlot(nDOFs,tBathe,1,cl1_POMEGRANATE,...
                   '-',cl1_POMEGRANATE,'s',0.5);
hMDED=customPlot(nDOFs,tMDED,1,cl3_DODGERBLUE,...
                   '--',cl3_DODGERBLUE,'o',1.0);

ylabel('CPU Time $$(s)$$','interpreter', 'latex')               
xlabel('Number of DOFs')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',12)
set(gca,'box','on')
set(gcf, 'Units','centimeters', 'Position',[5 5 16 8])
axis([0 5800 0 650])
grid on              
grid minor

legend([hBathe,hMDED],...     
     'Bathe',...
     'MDED','Location','northwest','Orientation','vertical','Position',[395 168 0.1 0.2])
 