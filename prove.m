clear;
clc;
%% initialization
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\FF100weightMV.mat'); % dim=100  
t = 100; % dimension
%*******************************************weight>0 || <0
% selected randomly from weightMV
matSize = size(weightMV);
randNum1 = unidrnd(matSize(1,2));
randNum2 = unidrnd(matSize(1,2));
w1 = weightMV(1:end,randNum1); % weight No.1
w2 = weightMV(1:end,randNum2); % weight No.2
%*******************************************weight>0 || <0

%*******************************************weight>0
% w1 = rand(t,1); % weight No.1
% w2 = rand(t,1); % weight No.2
% s1 = sum(w1); % sum of weight No.1
% s2 = sum(w2); % sum of weight No.2
% w1 = w1/s1; % scaling 0_1
% w2 = w2/s2; % scaling 0_1
%*******************************************weight>0
it = 100; % iteration time
r1 = ones(1,it); % result No.1 
                % define dim
r2 = ones(1,it); % result No.2
                % define dim
r = ones(1,it); % blending result
               % define dim
for i = 1:it
    rawData = rand(1, t+1);
    data = rawData(1,2:t+1)./rawData(1,1:t);
    for u = 0:0.01:1% blending coefficient
        w = u*w1 + (1-u)*w2;% blending weight
        r1(1,i) = data*w1;
        r2(1,i) = data*w2;
        r(1,i) = data*w;
       %% graph
        plot(r(1,it),'y');
        hold on;
        plot(r2(1,it),'r');
        hold on;
        plot(r1(1,it),'k');
        title('Simulation Result','fontsize',12);
        xlabel('Dimension','fontsize',12);
        ylabel('result','fontsize',12);
        legend('r','r2','r1');
    end
end
