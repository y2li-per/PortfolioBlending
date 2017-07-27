clear;
clc;
%% initialization
%load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\equity200_input.mat'); %weekly update
%load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ETF139_input.mat'); %weekly update
%load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff25_input.mat'); %monthly update
%load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff48_input.mat'); %monthly update
load('C:\Users\think\Documents\MATLAB\PortfolioBlending\data\ff100_input.mat'); %monthly update
%price = 1 + rand(100,701);
%gr = price(1:end,2:end)./price(1:end,1:end-1); % gross return 
                                               % first t for training
                                               % last m for testing
gr = portfolios.';
n = Nportfolios; %number of risky assets
m = Nmonths; %time period
             %498 months for FF
             %252 weeks for realworld
if (m == 498)
    t = 120; %training time period
             %120 months for FF
             %200 weeks for realworld
else
    t = 200;
end        
m = m - t; %investment period
%I = ones(n,1); % column vector all 1
weightEW = ones(n,m)/n; %weight vector of Equally-Weighted
weightVW = ones(n,m)/n; %weighy of Valued-Weighted(define dim)
weightMV = ones(n,m)/n; %weighy of Minimun-Variance(define dim)
weightMVVW = ones(n,m)/n; %weighy of MV-VW(define dim)
weightMVEW = ones(n,m)/n; %weighy of MV-EW(define dim)
weightEWVW = ones(n,m)/n; %weighy of EW-VW(define dim)
%% main algorithm
for k = 1:m
    %% computation of weightVW
    if(k == 1)
        weightVW(1:end,1) = gr(1:end,t)/sum(gr(1:end,t));
    else
        weightVW(1:end,k) = (weightVW(1:end,k-1).*gr(1:end,k+t-1))/(weightVW(1:end,k-1).'*gr(1:end,k+t-1));
        %sum(weightVW) == 1
    end
    %% computation of weightMV
    c = cov(gr(1:end,k:k+t-1).');
    weightMV(1:end,k) = (sum(inv(c))/sum(sum(inv(c)))).'; %sum(weightMV) == 1
    %% computation ofweightMVVW
    %if k == 1 || mod(k,10) == 0
        %if (k == 1)
            aMVVW = 1;%initialize Beta distribution Beta(a,b)
            bMVVW = 1;
        %end
        for j = 1:t
            aveU = aMVVW / (aMVVW+bMVVW);%average blending coefficient
            aveWeightMVVW = aveU*weightMV + (1-aveU)*weightVW; %average weighy of MV-VW
            %aveWeightMVEW = aveU*weightMV + (1-aveU)*weightEW; %average weighy of MV-EW
            aveRetMVVW = sum(aveWeightMVVW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-VW at time k
            %aveRetMVEW = sum(aveWeightMVEW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-EW at time k
            while (1)
                samU = betarnd(aMVVW,bMVVW); %sample blending coefficient
                if (samU ~= aveU) %resample when sample is equal to average
                    break;
                end
            end
            samWeightMVVW = samU*weightMV + (1-samU)*weightVW; %sample weighy of MV-VW
            %samWeightMVEW = samU*weightMV + (1-samU)*weightEW; %sample weighy of MV-EW
            samRetMVVW = sum(samWeightMVVW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-VW at time k
            %samRetMVEW = sum(samWeightMVEW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-EW at time k
            if (aveRetMVVW == samRetMVVW)
                %if (aveRetMVEW == samRetMVEW)
                flag = 3;
            elseif (((aveRetMVVW>samRetMVVW)&&(aveU>samU)) || ((aveRetMVVW<samRetMVVW)&&(aveU<samU)))
                %elseif (((aveRetMVEW>samRetMVEW)&&(aveU>samU)) || ((aveRetMVEW<samRetMVEW)&&(aveU<samU)))
                flag = 1;
            else
                flag = 2;
            end
            if (flag == 1) %success
                aMVVW = aMVVW + 1;
            elseif(flag == 2)%fail
                bMVVW = bMVVW + 1;
            end
        end
        u = (1+aMVVW)/(2+aMVVW+bMVVW); %decide final coefficient
        weightMVVW(1:end,k) = u*weightMV(1:end,k) + (1-u)*weightVW(1:end,k); %weighy of MV-VW(define dim)
        %weightMVEW(1:end,k) = u*weightMV(1:end,k) + (1-u)*weightEW(1:end,k); %weighy of MV-EW(define dim)
        %% computation of weightMVEW
        %if (k == 1)
            aMVEW = 1;%initialize Beta distribution Beta(a,b)
            bMVEW = 1;
        %end
        for j = 1:t
            aveU = aMVEW / (aMVEW+bMVEW);%average blending coefficient
            aveWeightMVEW = aveU*weightMV + (1-aveU)*weightEW; %average weighy of MV-EW
            aveRetMVEW = sum(aveWeightMVEW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-EW at time k
            while (1)
                samU = betarnd(aMVEW,bMVEW); %sample blending coefficient
                if (samU ~= aveU) %resample when sample is equal to average
                    break;
                end
            end
            samWeightMVEW = samU*weightMV + (1-samU)*weightEW; %sample weighy of MV-EW
            samRetMVEW = sum(samWeightMVEW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-EW at time k
            if (aveRetMVEW == samRetMVEW)
                flag = 3;
            elseif (((aveRetMVEW>samRetMVEW)&&(aveU>samU)) || ((aveRetMVEW<samRetMVEW)&&(aveU<samU)))
                flag = 1;
            else
                flag = 2;
            end
            if (flag == 1) %success
                aMVEW = aMVEW + 1;
            elseif(flag == 2)%fail
                bMVEW = bMVEW + 1;
            end
        end
        u = (1+aMVEW)/(2+aMVEW+bMVEW); %decide final coefficient
        weightMVEW(1:end,k) = u*weightMV(1:end,k) + (1-u)*weightEW(1:end,k); %weighy of MV-EW(define dim)
        %% computation of weightEWVW
        %if (k == 1)
            aEWVW = 1;%initialize Beta distribution Beta(a,b)
            bEWVW = 1;
        %end
        for j = 1:t
            aveU = aEWVW / (aEWVW+bEWVW);%average blending coefficient
            aveWeightEWVW = aveU*weightEW + (1-aveU)*weightVW; %average weighy of EW-VW
            aveRetEWVW = sum(aveWeightEWVW(1:end,k).*gr(1:end,j+k-1)); %return of ave EW-VW at time k
            while (1)
                samU = betarnd(aEWVW,bEWVW); %sample blending coefficient
                if (samU ~= aveU) %resample when sample is equal to average
                    break;
                end
            end
            samWeightEWVW = samU*weightEW + (1-samU)*weightVW; %sample weighy of MV-EW
            samRetEWVW = sum(samWeightEWVW(1:end,k).*gr(1:end,j+k-1)); %return of ave MV-EW at time k
            if (aveRetEWVW == samRetEWVW)
                flag = 3;
            elseif (((aveRetEWVW>samRetEWVW)&&(aveU>samU)) || ((aveRetEWVW<samRetEWVW)&&(aveU<samU)))
                flag = 1;
            else
                flag = 2;
            end
            if (flag == 1) %success
                aEWVW = aEWVW + 1;
            elseif(flag == 2)%fail
                bEWVW = bEWVW + 1;
            end
        end
        u = (1+aEWVW)/(2+aEWVW+bEWVW); %decide final coefficient
        weightEWVW(1:end,k) = u*weightEW(1:end,k) + (1-u)*weightVW(1:end,k); %weighy of MV-EW(define dim)
%     else
%         weightMVEW(1:end,k) = weightMVEW(1:end,k-1);
%         weightMVVW(1:end,k) = weightMVVW(1:end,k-1);
%         weightEWVW(1:end,k) = weightEWVW(1:end,k-1);
%     end
end

%% computation of average change rate
chaRaEW = weightEW(1:end,2:end); % define dim
chaRaVW = weightVW(1:end,2:end); % define dim
chaRaMV = weightMV(1:end,2:end); % define dim
chaRaMVEW = weightMVEW(1:end,2:end); % define dim
chaRaMVVW = weightMVVW(1:end,2:end); % define dim
chaRaEWVW = weightEWVW(1:end,2:end); % define dim
for i = 1:m-1
    chaRaEW(1:end,i) = abs((chaRaEW(1:end,i)-weightEW(1:end,i))./weightEW(1:end,i));
    chaRaVW(1:end,i) = abs((chaRaVW(1:end,i)-weightVW(1:end,i))./weightVW(1:end,i));
    chaRaMV(1:end,i) = abs((chaRaMV(1:end,i)-weightMV(1:end,i))./weightMV(1:end,i));
    chaRaMVEW(1:end,i) = abs((chaRaMVEW(1:end,i)-weightMVEW(1:end,i))./weightMVEW(1:end,i));
    chaRaMVVW(1:end,i) = abs((chaRaMVVW(1:end,i)-weightMVVW(1:end,i))./weightMVVW(1:end,i));
    chaRaEWVW(1:end,i) = abs((chaRaEWVW(1:end,i)-weightEWVW(1:end,i))./weightEWVW(1:end,i));
end
aveChaRaEW = mean(mean(chaRaEW));
aveChaRaVW = mean(mean(chaRaVW));
aveChaRaMV = mean(mean(chaRaMV));
aveChaRaMVEW = mean(mean(chaRaMVEW));
aveChaRaMVVW = mean(mean(chaRaMVVW));
aveChaRaEWVW = mean(mean(chaRaEWVW));
%% 
%attach to ThompsonSampling.m when using
%or run ThompsonSampling.m first
%% computation of return at each period
%**********************************return matrics
allRetEW = weightEW.*gr(1:end,t+1:end);
allRetVW = weightVW.*gr(1:end,t+1:end);
allRetMV = weightMV.*gr(1:end,t+1:end);
allRetMVVW = weightMVVW.*gr(1:end,t+1:end);
allRetMVEW = weightMVEW.*gr(1:end,t+1:end);
allRetEWVW = weightEWVW.*gr(1:end,t+1:end);
%**********************************return matrics

%**********************************return vector
retEW = sum(allRetEW);
retVW = sum(allRetVW);
retMV = sum(allRetMV);
retMVVW = sum(allRetMVVW);
retMVEW = sum(allRetMVEW);
retEWVW = sum(allRetEWVW);
%**********************************return vector

%**********************************net return vector
netRetEW = retEW - 1;
netRetVW = retVW - 1;
netRetMV = retMV - 1;
netRetMVVW = retMVVW - 1;
netRetMVEW = retMVEW - 1;
netRetEWVW = retEWVW - 1;
%**********************************net return vector
%% computation of Sharp Ratio
costR = 1/50; %cost factor
%**********************************EW
tmpWeiEW = weightEW; %define dim
tmpWeiEW(1:end,2:end) = weightEW(1:end,1:end-1);
afNetRetEW = netRetEW.*(1-costR*sum(abs(weightEW-tmpWeiEW))); %after-cost net return
aveAfNetRetEW = mean(afNetRetEW); %average after-cost net return
sdEW = std(afNetRetEW,1); %standard deviation
srEW = aveAfNetRetEW/sdEW; %sharp ratio
%**********************************EW

%**********************************VW
tmpWeiVW = weightVW; %define dim
tmpWeiVW(1:end,2:end) = weightVW(1:end,1:end-1);
afNetRetVW = netRetVW.*(1-costR*sum(abs(weightVW-tmpWeiVW))); %after-cost net return
aveAfNetRetVW = mean(afNetRetVW); %average after-cost net return
sdVW = std(afNetRetVW,1); %standard deviation
srVW = aveAfNetRetVW/sdVW; %sharp ratio
%**********************************VW

%**********************************MV
tmpWeiMV = weightMV; %define dim
tmpWeiMV(1:end,2:end) = weightMV(1:end,1:end-1);
afNetRetMV = netRetMV.*(1-costR*sum(abs(weightMV-tmpWeiMV))); %after-cost net return
aveAfNetRetMV = mean(afNetRetMV); %average after-cost net return
sdMV = std(afNetRetMV,1); %standard deviation
srMV = aveAfNetRetMV/sdMV; %sharp ratio
%**********************************MV

%**********************************MVEW
tmpWeiMVEW = weightMVEW; %define dim
tmpWeiMVEW(1:end,2:end) = weightMVEW(1:end,1:end-1);
afNetRetMVEW = netRetMVEW.*(1-costR*sum(abs(weightMVEW-tmpWeiMVEW))); %after-cost net return
aveAfNetRetMVEW = mean(afNetRetMVEW); %average after-cost net return
sdMVEW = std(afNetRetMVEW,1); %standard deviation
srMVEW = aveAfNetRetMVEW/sdMVEW; %sharp ratio
%**********************************MVEW

%**********************************MVVW
tmpWeiMVVW = weightMVVW; %define dim
tmpWeiMVVW(1:end,2:end) = weightMVVW(1:end,1:end-1);
afNetRetMVVW = netRetMVVW.*(1-costR*sum(abs(weightMVVW-tmpWeiMVVW))); %after-cost net return
aveAfNetRetMVVW = mean(afNetRetMVVW); %average after-cost net return
sdMVVW = std(afNetRetMVVW,1); %standard deviation
srMVVW = aveAfNetRetMVVW/sdMVVW; %sharp ratio
%**********************************MVVW

%**********************************EWVW
tmpWeiEWVW = weightEWVW; %define dim
tmpWeiEWVW(1:end,2:end) = weightEWVW(1:end,1:end-1);
afNetRetEWVW = netRetEWVW.*(1-costR*sum(abs(weightEWVW-tmpWeiEWVW))); %after-cost net return
aveAfNetRetEWVW = mean(afNetRetEWVW); %average after-cost net return
sdEWVW = std(afNetRetEWVW,1); %standard deviation
srEWVW = aveAfNetRetEWVW/sdEWVW; %sharp ratio
%**********************************EWVW
%% computation of volatility
if (m == 498-t)
    H = 12; %total number of rebalancing times each year
             %12 months/year for FF
             %52 weeks/year for realworld
else
    H = 52;
end 

vEW = sqrt(H) * sdEW; %volatility
vVW = sqrt(H) * sdVW; %volatility
vMV = sqrt(H) * sdMV; %volatility
vMVEW = sqrt(H) * sdMVEW; %volatility
vMVVW = sqrt(H) * sdMVVW; %volatility
vEWVW = sqrt(H) * sdEWVW; %volatility
%% computation of Maximum Drawdown
%**********************************EW
cuNetRetEW = sum((allRetEW).*(1-costR*abs(weightEW-tmpWeiEW))); %after-cost net return
for i = 2:m
    cuNetRetEW(1,i) = cuNetRetEW(1,i)*cuNetRetEW(1,i-1);
end
detEW = cuNetRetEW; %define dim
for i = 1:m
    detEW(1,i) = (max(cuNetRetEW(1,1:i)) - cuNetRetEW(1,i))/max(cuNetRetEW(1,1:i));
end
mddEW = max(detEW(1,1:end));
%**********************************EW

%**********************************VW
cuNetRetVW = sum((allRetVW).*(1-costR*abs(weightVW-tmpWeiVW))); %after-cost net return
for i = 2:m
    cuNetRetVW(1,i) = cuNetRetVW(1,i)*cuNetRetVW(1,i-1);
end
detVW = cuNetRetVW; %define dim
for i = 1:m   
    detVW(1,i) = (max(cuNetRetVW(1,1:i)) - cuNetRetVW(1,i))/max(cuNetRetVW(1,1:i));
end
mddVW = max(detVW(1,1:end));
%**********************************VW

%**********************************MV
cuNetRetMV = sum((allRetMV).*(1-costR*abs(weightMV-tmpWeiMV))); %after-cost net return
for i = 2:m
    cuNetRetMV(1,i) = cuNetRetMV(1,i)*cuNetRetMV(1,i-1);
end
detMV = cuNetRetMV; %define dim
for i = 1:m
    detMV(1,i) = (max(cuNetRetMV(1,1:i)) - cuNetRetMV(1,i))/max(cuNetRetMV(1,1:i));
end
mddMV = max(detMV(1,1:end));
%**********************************MV

%**********************************MVEW
cuNetRetMVEW = sum((allRetMVEW).*(1-costR*abs(weightMVEW-tmpWeiMVEW))); %after-cost net return
for i = 2:m
    cuNetRetMVEW(1,i) = cuNetRetMVEW(1,i)*cuNetRetMVEW(1,i-1);
end
detMVEW = cuNetRetMVEW; %define dim
for i = 1:m
    detMVEW(1,i) = (max(cuNetRetMVEW(1,1:i)) - cuNetRetMVEW(1,i))/max(cuNetRetMVEW(1,1:i));
end
mddMVEW = max(detMVEW(1,1:end));
%**********************************MVEW

%**********************************MVVW
cuNetRetMVVW = sum((allRetMVVW).*(1-costR*abs(weightMVVW-tmpWeiMVVW))); %after-cost net return
for i = 2:m
    cuNetRetMVVW(1,i) = cuNetRetMVVW(1,i)*cuNetRetMVVW(1,i-1);
end
detMVVW = cuNetRetMVVW; %define dim
for i = 1:m
    detMVVW(1,i) = (max(cuNetRetMVVW(1,1:i)) - cuNetRetMVVW(1,i))/max(cuNetRetMVVW(1,1:i));
end
mddMVVW = max(detMVVW(1,1:end));
%**********************************MVVW

%**********************************EWVW
cuNetRetEWVW = sum((allRetEWVW).*(1-costR*abs(weightEWVW-tmpWeiEWVW))); %after-cost net return
for i = 2:m
    cuNetRetEWVW(1,i) = cuNetRetEWVW(1,i)*cuNetRetEWVW(1,i-1);
end
detEWVW = cuNetRetEWVW; %define dim
for i = 1:m
    detEWVW(1,i) = (max(cuNetRetEWVW(1,1:i)) - cuNetRetEWVW(1,i))/max(cuNetRetEWVW(1,1:i));
end
mddEWVW = max(detEWVW(1,1:end));
%**********************************EWVW
%% computation of cumulative return
cuRetEW = retEW; %define dim
cuRetVW = retVW; %define dim
cuRetMV = retMV; %define dim
cuRetMVVW = retMVVW; %define dim
cuRetMVEW = retMVEW; %define dim
cuRetEWVW = retEWVW; %define dim
for i = 2:m
    cuRetEW(1,i) = retEW(1,i)*cuRetEW(1,i-1);
    cuRetVW(1,i) = retVW(1,i)*cuRetVW(1,i-1);
    cuRetMV(1,i) = retMV(1,i)*cuRetMV(1,i-1);
    cuRetMVVW(1,i) = retMVVW(1,i)*cuRetMVVW(1,i-1);
    cuRetMVEW(1,i) = retMVEW(1,i)*cuRetMVEW(1,i-1);
    cuRetEWVW(1,i) = retEWVW(1,i)*cuRetEWVW(1,i-1);
end
%% graph of EW
plot(cuRetEW,'g');
title('EW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Cumulative Wealth','fontsize',12);
hold on;
%% graph of VW
plot(cuRetVW,'r');
title('VW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Cumulative Wealth','fontsize',12);
hold on;
%% graph of MV
plot(cuRetMV,'c');
title('MV','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Cumulative Wealth','fontsize',12);
hold on;
%% graph of MVVW
plot(cuRetMVVW,'b');
title('MV-VW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
hold on;
%% graph of MVEW
plot(cuRetMVEW,'m');
title('MV-EW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Cumulative Wealth','fontsize',12);
hold on;
%% graph of EWVW
plot(cuRetEWVW,'y');
title('EW-VW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Cumulative Wealth','fontsize',12);
%% 
title('Blending Portfolio Performance by EQ181','fontsize',20);
legend('Equally-Weight','Value-Weighted','Minimum-Variance','MV-VW','MV-EW','EW-VW');



