clear;
clc;

ETF1_main = readmatrix('IVE.xlsx');
ETF2_main = readmatrix('IAU.csv');
ETF3_main = readmatrix('LQD.csv');
SP = readmatrix('SP500.csv');


assetRetns(:,1) = log(ETF1_main(2:2518,6)./ETF1_main(1:2517,6));
assetRetns(:,2) = log(ETF2_main(2:2518,6)./ETF2_main(1:2517,6));
assetRetns(:,3) = log(ETF3_main(2:2518,6)./ETF3_main(1:2517,6));
benchRetn = SP(2:2518,6);


%% MV Frontier
Sigma = cov(assetRetns);
port = Portfolio('NumAssets', numAssets, 'lb', 0, 'budget', 1, 'Name', 'Mean Variance');
port = setAssetMoments(port, mean(assetRetns), Sigma);
wts = estimateMaxSharpeRatio(port);
[risk, ret] = estimatePortMoments(port, wts);
plotFrontier(port, 20);
hold on
plot(risk,ret,'*r');
hold off

%% Black Litterman
assetNames = ["IVE", "IAU", "LQD"];
numAssets = size(assetRetns, 2);
v = 2;
P = zeros(v, numAssets);
q = zeros(v, 1);
Omega = zeros(v);

% View 1
P(1, assetNames=="IAU") = 1; 
q(1) = 0.25;
Omega(1, 1) = 1e-5;

% View 2
P(2, assetNames=="IVE") = 1; 
q(2) = 0.09;
Omega(2, 2) = 1e-2;

viewTable = array2table([P q diag(Omega)], 'VariableNames', [assetNames "View_Return" "View_Uncertainty"]) 

bizyear2bizday = 1/252;
q = q*bizyear2bizday; 
Omega = Omega*bizyear2bizday;

tau = 1/size(assetRetns,1);
C = tau*Sigma;

[wtsMarket, PI] = findMarketPortfolioAndImpliedReturn(assetRetns, benchRetn);
mu_bl = (P'*(Omega\P) + inv(C)) \ ( C\PI + P'*(Omega\q));
cov_mu = inv(P'*(Omega\P) + inv(C));
table(assetNames', PI*252, mu_bl*252, 'VariableNames', ["Asset_Name", ...
    "Prior_Belief_of_Expected_Return", "Black_Litterman_Blended_Expected_Return"])


portBL = Portfolio('NumAssets', numAssets, 'lb', 0, 'budget', 1, 'Name', 'Mean Variance with Black-Litterman');
portBL = setAssetMoments(portBL, mu_bl, Sigma + cov_mu);  
wtsBL = estimateMaxSharpeRatio(portBL);
[risk, ret] = estimatePortMoments(portBL, wtsBL);
plotFrontier(portBL, 20);
hold on
plot(risk,ret,'*r');

%% GARCH 
[portG,weightGARCH] = getGarchData();


%% Construct pie chart of weight distribution
ax1 = subplot(2,2,1);
idx = wts>0.001;
pie(ax1, wts(idx), assetNames(idx));
title(ax1, port.Name ,'Position', [-0.05, 1.5, 0]);

ax2 = subplot(2,2,2);
idx_BL = wtsBL>0.001;
pie(ax2, wtsBL(idx_BL), assetNames(idx_BL));
title(ax2, portBL.Name ,'Position', [-0.05, 1.5, 0]);

ax3 = subplot(2,2,3);
idx_G = weightGARCH>0.001;
pie(ax3, weightGARCH(idx_G), assetNames(idx_G));
title(ax3, portG.Name ,'Position', [-0.05, 1.5, 0]);


function [wtsMarket, PI] = findMarketPortfolioAndImpliedReturn(assetRetn, benchRetn)
    Sigma = cov(assetRetn);
    numAssets = size(assetRetn,2);
    LB = zeros(1,numAssets);
    Aeq = ones(1,numAssets);
    Beq = 1;
    opts = optimoptions('lsqlin','Algorithm','interior-point', 'Display',"off");
    wtsMarket = lsqlin(assetRetn, benchRetn, [], [], Aeq, Beq, LB, [], [], opts);
    shpr = mean(benchRetn)/std(benchRetn);
    delta = shpr/sqrt(wtsMarket'*Sigma*wtsMarket); 
    PI = delta*Sigma*wtsMarket;
end


function [port,wtsGARCH] = getGarchData()
ETF1_main = readmatrix('IVE.xlsx');
ETF2_main = readmatrix('IAU.csv');
ETF3_main = readmatrix('LQD.csv');

ETF_test(:,1) = ETF1_main(:,6);
ETF_test(:,2) = ETF2_main(:,6);
ETF_test(:,3) = ETF3_main(:,6);


ETF(:,1) = log(ETF1_main(2:2518,6)./ETF1_main(1:2517,6));
ETF(:,2) = log(ETF2_main(2:2518,6)./ETF2_main(1:2517,6));
ETF(:,3) = log(ETF3_main(2:2518,6)./ETF3_main(1:2517,6));



EstMdl11 = estimate(garch(1,1),ETF(:,1));
EstMdl12 = estimate(garch(1,2),ETF(:,1));
EstMdl13 = estimate(garch(1,3),ETF(:,1));
EstMdl21 = estimate(garch(2,1),ETF(:,1));
EstMdl22 = estimate(garch(2,2),ETF(:,1));
EstMdl23 = estimate(garch(2,3),ETF(:,1));
EstMdl33 = estimate(garch(3,3),ETF(:,1));
EstMdl31 = estimate(garch(3,1),ETF(:,1));
EstMdl32 = estimate(garch(3,2),ETF(:,1));


results= summarize(EstMdl11);
mdleval(1,1) = results.AIC;
mdleval(1,2) = results.BIC;
results= summarize(EstMdl12);
mdleval(2,1) = results.AIC;
mdleval(2,2) = results.BIC;
results= summarize(EstMdl13);
mdleval(3,1) = results.AIC;
mdleval(3,2) = results.BIC;
results= summarize(EstMdl21);
mdleval(4,1) = results.AIC;
mdleval(4,2) = results.BIC;
results= summarize(EstMdl22);
mdleval(5,1) = results.AIC;
mdleval(5,2) = results.BIC;
results= summarize(EstMdl23);
mdleval(6,1) = results.AIC;
mdleval(6,2) = results.BIC;
results= summarize(EstMdl31);
mdleval(7,1) = results.AIC;
mdleval(7,2) = results.BIC;
results= summarize(EstMdl32);
mdleval(8,1) = results.AIC;
mdleval(8,2) = results.BIC;
results= summarize(EstMdl33);
mdleval(9,1) = results.AIC;
mdleval(9,2) = results.BIC;
[minVal, index] = min(mdleval);

garchPred(:,1) = infer(EstMdl22,ETF(:,1));

EstMdl11 = estimate(garch(1,1),ETF(:,2));
EstMdl12 = estimate(garch(1,2),ETF(:,2));
EstMdl13 = estimate(garch(1,3),ETF(:,2));
EstMdl21 = estimate(garch(2,1),ETF(:,2));
EstMdl22 = estimate(garch(2,2),ETF(:,2));
EstMdl23 = estimate(garch(2,3),ETF(:,2));
EstMdl33 = estimate(garch(3,3),ETF(:,2));
EstMdl31 = estimate(garch(3,1),ETF(:,2));
EstMdl32 = estimate(garch(3,2),ETF(:,2));

results= summarize(EstMdl11);
mdleval(1,1) = results.AIC;
mdleval(1,2) = results.BIC;
results= summarize(EstMdl12);
mdleval(2,1) = results.AIC;
mdleval(2,2) = results.BIC;
results= summarize(EstMdl13);
mdleval(3,1) = results.AIC;
mdleval(3,2) = results.BIC;
results= summarize(EstMdl21);
mdleval(4,1) = results.AIC;
mdleval(4,2) = results.BIC;
results= summarize(EstMdl22);
mdleval(5,1) = results.AIC;
mdleval(5,2) = results.BIC;
results= summarize(EstMdl23);
mdleval(6,1) = results.AIC;
mdleval(6,2) = results.BIC;
results= summarize(EstMdl31);
mdleval(7,1) = results.AIC;
mdleval(7,2) = results.BIC;
results= summarize(EstMdl32);
mdleval(8,1) = results.AIC;
mdleval(8,2) = results.BIC;
results= summarize(EstMdl33);
mdleval(9,1) = results.AIC;
mdleval(9,2) = results.BIC;
[minVal, index] = min(mdleval);

garchPred(:,2) = infer(EstMdl11,ETF(:,2));

EstMdl11 = estimate(garch(1,1),ETF(:,3));
EstMdl12 = estimate(garch(1,2),ETF(:,3));
EstMdl13 = estimate(garch(1,3),ETF(:,3));
EstMdl21 = estimate(garch(2,1),ETF(:,3));
EstMdl22 = estimate(garch(2,2),ETF(:,3));
EstMdl23 = estimate(garch(2,3),ETF(:,3));
EstMdl33 = estimate(garch(3,3),ETF(:,3));
EstMdl31 = estimate(garch(3,1),ETF(:,3));
EstMdl32 = estimate(garch(3,2),ETF(:,3));

results= summarize(EstMdl11);
mdleval(1,1) = results.AIC;
mdleval(1,2) = results.BIC;
results= summarize(EstMdl12);
mdleval(2,1) = results.AIC;
mdleval(2,2) = results.BIC;
results= summarize(EstMdl13);
mdleval(3,1) = results.AIC;
mdleval(3,2) = results.BIC;
results= summarize(EstMdl21);
mdleval(4,1) = results.AIC;
mdleval(4,2) = results.BIC;
results= summarize(EstMdl22);
mdleval(5,1) = results.AIC;
mdleval(5,2) = results.BIC;
results= summarize(EstMdl23);
mdleval(6,1) = results.AIC;
mdleval(6,2) = results.BIC;
results= summarize(EstMdl31);
mdleval(7,1) = results.AIC;
mdleval(7,2) = results.BIC;
results= summarize(EstMdl32);
mdleval(8,1) = results.AIC;
mdleval(8,2) = results.BIC;
results= summarize(EstMdl33);
mdleval(9,1) = results.AIC;
mdleval(9,2) = results.BIC;
[minVal, index] = min(mdleval);


garchPred(:,3) = infer(EstMdl31,ETF(:,3));

garchMean = mean(garchPred);
garchCov = cov(garchPred);

port = Portfolio('NumAssets', 3, 'lb', 0, 'budget', 1, 'Name', 'Garch Model');
port = setAssetMoments(port, garchMean, garchCov);
plotFrontier(port, 20);
wtsGARCH = estimateMaxSharpeRatio(port);
[risk, ret] = estimatePortMoments(port, wtsGARCH);

end
