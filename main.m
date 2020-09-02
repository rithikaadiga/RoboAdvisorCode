clear;
clc;

%% Load Data
ETF1_main = readmatrix('IVE.xlsx');
ETF2_main = readmatrix('IAU.csv');
ETF3_main = readmatrix('LQD.csv');
ETF(:,1) = ETF1_main(2:2518,1);ETF(:,1)
dates = x2mdate(ETF(:,1),0,'datetime')

ETF_stats = [ETF1_main(:,6) ETF2_main(:,6) ETF3_main(:,6)];

% Calculating log normal returns
ETF(:,1) = log(ETF1_main(2:2518,6)./ETF1_main(1:2517,6));
ETF(:,2) = log(ETF2_main(2:2518,6)./ETF2_main(1:2517,6));
ETF(:,3) = log(ETF3_main(2:2518,6)./ETF3_main(1:2517,6));


%% Descriptive Stats
retMeans = mean(ETF_stats);
retStd = std(ETF_stats);
retVarCov = cov(ETF_stats);
retMedian = median(ETF_stats);
retKurt = kurtosis(ETF_stats);
retSkewness = skewness(ETF_stats);
retRange = range(ETF_stats)
retCorr = corrcoef(ETF_stats);
rfAvg = mean(RF)/100;

%% Histogram of Normal Returns
histfit(ETF1_main(:,6))
title('Histogram with distribution fit for IVE');

histfit(ETF2_main(:,6))
title('Histogram with distribution fit for IAU');

histfit(ETF3_main(:,6))
title('Histogram with distribution fit for LQD');


%% JB Test : When jbstat > critical, the null hypothesis that X is normally distributed can be rejected
ETF_test(:,1) = ETF1_main(:,6);
ETF_test(:,2) = ETF2_main(:,6);
ETF_test(:,3) = ETF3_main(:,6);

[h,p,jbstat,critical] = jbtest(ETF_test(:,1),0.05);
[h,p,jbstat,critical] = jbtest(ETF_test(:,2),0.05);
[h,p,jbstat,critical] = jbtest(ETF_test(:,3),0.05);

%% Ljung-Box Test
returns = price2ret(ETF_test(:,1));
res = returns - mean(returns);
[h1,pvalue] = lbqtest(res);
% reject the null hypothesis that the residuals are not autocorrelated.
returns = price2ret(ETF_test(:,2));
res = returns - mean(returns);
[h2,pvalue] = lbqtest(res);
% may not be able to reject here
returns = price2ret(ETF_test(:,3));
res = returns - mean(returns);
[h3,pvalue] = lbqtest(res);
% reject the null hypothesis that the residuals are not autocorrelated.

%% Augmented Dickey Fuller Test
hdf1 = adftest(ETF_test(:,1),'lags',0:3)
hdf2 = adftest(ETF_test(:,2),'lags',0:3)
hdf3 = adftest(ETF_test(:,3),'lags',0:3)
%  all three tests fail to reject the null hypothesis of a unit root

%% Time vs Log Returns
hold on;
plot(dates,ETF(:,1))
title('Log Normal Returns of iShares S&P 500 Value Index (IVE)');
xlabel('Time');
ylabel('Log Returns');
hold off;

hold on;
plot(dates,ETF(:,2));
title('Log Normal Returns of iShares Gold Trust (IAU)');
xlabel('Time');
ylabel('Log Returns');
hold off;

hold on;
plot(dates,ETF(:,3))
title('Log Normal Returns of iShares iBoxx $ Investment Grade Corporate Bond ETF (LDQ)');
xlabel('Time');
ylabel('Log Returns');
hold off;

%% Histogram of Log-Normal Returns
histfit(ETF(:,1))
title('Histogram with distribution fit for IVE');

histfit(ETF(:,2))
title('Histogram with distribution fit for IAU');

histfit(ETF(:,3))
title('Histogram with distribution fit for LQD');