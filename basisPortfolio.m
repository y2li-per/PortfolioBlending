%attach to ThompsonSampling.m when using
%or run ThompsonSampling.m first
%************************************************computation of return at each period
retEW = sum(weightEW.*gr(1:end,t+1:end));
retVW = sum(weightVW.*gr(1:end,t+1:end));
retMV = sum(weightMV.*gr(1:end,t+1:end));
retMVVW = sum(weightMVVW.*gr(1:end,t+1:end));
retMVEW = sum(weightMVEW.*gr(1:end,t+1:end));
%************************************************computation of return at each period

%************************************************computation of acumulative return
for i = 2:m
    retEW(1,i) = retEW(1,i)*retEW(1,i-1);
    retVW(1,i) = retVW(1,i)*retVW(1,i-1);
    retMV(1,i) = retMV(1,i)*retMV(1,i-1);
    retMVVW(1,i) = retMVVW(1,i)*retMVVW(1,i-1);
    retMVEW(1,i) = retMVEW(1,i)*retMVEW(1,i-1);
end
%************************************************computation of acumulative return

%************************************************graph of EW
subplot(2,3,1);
plot(retEW,'b');
title('EW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
hold on;
%************************************************graph of EW

%************************************************graph of VW
subplot(2,3,2);
plot(retVW,'r');
title('VW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
hold on;
%************************************************graph of VW

%************************************************graph of MV
subplot(2,3,3);
plot(retMV,'g');
title('MV','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
%************************************************graph of MV

%************************************************graph of MVVW
subplot(2,3,4);
plot(retMVVW,'b');
title('MV-VW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
%************************************************graph of MVVW

%************************************************graph of MVEW
subplot(2,3,6);
plot(retMV,'r');
title('MV-EW','fontsize',12);
xlabel('Investment Period','fontsize',12);
ylabel('Acumulative Wealth','fontsize',12);
%************************************************graph of MVEW

title('Basis Portfolio Performance by FF100','fontsize',20);
legend('Equally-Weight','Value-Weighted','Minimum-Variance')