clear;
clc;
%% loading
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff100_input.mat'); %monthly update
%% graph of price
%subplot(1,5,1);
%plot(portfolios_price,'r');
%title('Price','fontsize',12);
%xlabel('Time Period','fontsize',12);
%ylabel('Price','fontsize',12);
%% graph of return
subplot(1,5,5);
plot(sum(portfolios),'b');
title('Return of FF100','fontsize',12);
xlabel('Time Period','fontsize',12);
ylabel('Return','fontsize',12);
grid on;
%% loading
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\equity200_input.mat'); %weekly update
%% graph of return
subplot(1,5,1);
plot(sum(portfolios),'b');
title('Return of EQ181','fontsize',12);
xlabel('Time Period','fontsize',12);
ylabel('Return','fontsize',12);
grid on;
%% loading
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ETF139_input.mat'); %weekly update
%% graph of return
subplot(1,5,2);
plot(sum(portfolios),'b');
title('Return of ETF139','fontsize',12);
xlabel('Time Period','fontsize',12);
ylabel('Return','fontsize',12);
grid on;
%% loading
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff25_input.mat'); %monthly update
%% graph of return
subplot(1,5,3);
plot(sum(portfolios),'b');
title('Return of FF25','fontsize',12);
xlabel('Time Period','fontsize',12);
ylabel('Return','fontsize',12);
grid on;
%% loading
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff48_input.mat'); %monthly update
%% graph of return
subplot(1,5,4);
plot(sum(portfolios),'b');
title('Return of FF48','fontsize',12);
xlabel('Time Period','fontsize',12);
ylabel('Return','fontsize',12);
grid on;
