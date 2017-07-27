clear;
clc;
%% initialization
dim = 100;
up = ones(1,dim); % upper confidence
                  % define dim
low = ones(1,dim); % lower confidence
                   % define dim
mid = ones(1,dim); % blending
                   % define dim
%% sampling
for i = 2:dim
    up(1,i) = 0.9+rand(1,1);
    while(up(1,i) > 1.2)
        up(1,i) = 0.8+rand(1,1);
    end
    mid(1,i) = 0.8+rand(1,1);
    while(mid(1,i) > up(1,i))
        mid(1,i) = 0.8+rand(1,1);
    end
    low(1,i) = 0.8+rand(1,1);
    while(low(1,i) > mid(1,i))
        low(1,i) = 0.8+rand(1,1);
    end
end
%% computation of cumulative result
rup = up; % result of up
          % define dim
rlow = low; % result of low
            % define dim
rmid = mid; % result of middle
            % define dim
for i = 2:dim
    rup(1,i) = rup(1,i-1) * up(1,i);
    rlow(1,i) = rlow(1,i-1) * low(1,i);
    rmid(1,i) = rmid(1,i-1) * mid(1,i);
end
%% graph
%***************************************graph1
subplot(1,2,1);
plot(up,'k');
hold on;
plot(low,'r');
hold on;
plot(mid,'y');
title('Simulation Result','fontsize',12);
xlabel('Time','fontsize',12);
ylabel('Weight','fontsize',12);
legend('up','low','middle');
%***************************************graph1

%***************************************graph2
subplot(1,2,2);
plot(rup,'k');
hold on;
plot(rlow,'r');
hold on;
plot(rmid,'y');
title('Simulation Result','fontsize',12);
xlabel('Time','fontsize',12);
ylabel('Cumulative Result','fontsize',12);
legend('up','low','middle');
 %***************************************graph2