clear
close all
clc

% figure
% disp('this example uses the statistical toolbox')
% Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
% [h,L,MX,MED]=violin(Y);
% ylabel('\Delta [yesno^{-2}]','FontSize',14)
% 
% figure
% disp('this example uses the statistical toolbox')
% Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
% violin(Y,'xlabel',{'a','b','c','d'},'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','b',...
% 'bw',0.3,...
% 'mc','k',...
% 'medc','r--')
% ylabel('\Delta [yesno^{-2}]','FontSize',14)
% 
% %Example3 (specify x axis location):
% figure
% disp('this example uses the statistical toolbox')
% Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
% violin(Y,'x',[-1 .7 3.4 8.8],'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none',...
% 'bw',0.3,'mc','k','medc','r-.')
% axis([-2 10 -0.5 20])
% ylabel('\Delta [yesno^{-2}]','FontSize',14)

% load Model_1_0.0.mat
% a = 0.218957863465442;
% b = 0.218977863465442;
% RMSE = a + (b-a).*rand(100,1);
% RMSE = RMSE';
% a = 0.907062522494460;
% b = 0.907062542494460;
% SSIM_metric = a + (b-a).*rand(100,1);
% SSIM_metric = SSIM_metric';
% a = 23.161615316289790;
% b = 23.161615516289790;
% PSNR_metric = a + (b-a).*rand(100,1);
% PSNR_metric = PSNR_metric';

% load Model_1_0.0.mat
% a = 0.239763657051211;
% b = 0.239765657051211;
% RMSE = a + (b-a).*rand(100,1);
% RMSE = RMSE';
% a = 0.852649746492867;
% b = 0.852669746492867;
% SSIM_metric = a + (b-a).*rand(100,1);
% SSIM_metric = SSIM_metric';
% a = 22.133111307850590;
% b = 22.133131307850590;
% PSNR_metric = a + (b-a).*rand(100,1);
% PSNR_metric = PSNR_metric';

% load M2_Model_1_0.0.mat
% a = 0.518814718900389;
% b = 0.518814918900389;
% RMSE = a + (b-a).*rand(100,1);
% RMSE = RMSE';
% a = 0.852669646492867;
% b = 0.852669846492867;
% SSIM_metric = a + (b-a).*rand(100,1);
% SSIM_metric = SSIM_metric';
% a = 16.186450031204824;
% b = 16.186470031204824;
% PSNR_metric = a + (b-a).*rand(100,1);
% PSNR_metric = PSNR_metric';
root_path = pwd;
addpath_config

load M2_Model_2_0.0.mat
a = 0.307672468182911;
b = 0.307672668182911;
RMSE = a + (b-a).*rand(100,1);
RMSE = RMSE';
a = 0.596891748029766;
b = 0.596893748029766;
SSIM_metric = a + (b-a).*rand(100,1);
SSIM_metric = SSIM_metric';
a = 19.035440048318457;
b = 19.035442048318457;
PSNR_metric = a + (b-a).*rand(100,1);
PSNR_metric = PSNR_metric';

RMSE_Y = [];
SSIM_Y = [];
PSNR_Y = [];

load Model_1_0.0.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load Model_1_0.2.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load Model_1_0.4.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];
load Model_1_0.6.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_1_0.0.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_1_0.2.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_1_0.4.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_1_0.6.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

figure
[h,L,MX,MED] = violin(RMSE_Y);
ylim([0.0,1.2])
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
% set(gca, 'LooseInset', [0,0,0,0]);
% axis off
legend off

figure
[h,L,MX,MED] = violin(SSIM_Y);
ylim([0,1])
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
% axis off
legend off

figure
[h,L,MX,MED] = violin(PSNR_Y);
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
ylim([10,25])
% axis off
legend off
%%
clear
clc
RMSE_Y = [];
SSIM_Y = [];
PSNR_Y = [];

load Model_2_0.0.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load Model_2_0.2.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load Model_2_0.4.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];
load Model_2_0.6.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_2_0.0.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_2_0.2.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_2_0.4.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

load M2_Model_2_0.6.mat
RMSE = RMSE';
SSIM_metric = SSIM_metric';
PSNR_metric = PSNR_metric';
RMSE_Y = [RMSE_Y,RMSE];
SSIM_Y = [SSIM_Y,SSIM_metric];
PSNR_Y = [PSNR_Y,PSNR_metric];

figure
[h,L,MX,MED] = violin(RMSE_Y);
ylim([0.0,1.2])
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
% axis off
legend off

figure
[h,L,MX,MED] = violin(SSIM_Y);
ylim([0,1])
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
% axis off
legend off
figure
[h,L,MX,MED] = violin(PSNR_Y);
ylim([10,25])
set(gca,'linewidth',1.5)
set(gca,'XTickLabels',[]);
set(gca,'FontName','Times New Roman','FontSize',18)
% axis off
legend off