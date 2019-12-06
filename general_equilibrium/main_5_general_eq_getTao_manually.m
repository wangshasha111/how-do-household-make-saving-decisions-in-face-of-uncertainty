%% Question 3.3.2
% Evaluate Universal Basic Income

%% Step 2
% Then evaluate Andrew Young＊s proposal of a UBI: Compare the stationary
% equilibrium with 竹 = 而 = 0 to a stationary equilibrium where 竹 = 0:2 and the
% tax rate 而 adjusts to clear the government budget. What happens to macroeconomic 
% aggregates (Y; K; C), to equilibrium prices (w; r) and the equilibrium
% distributions for earnings (1 ? 而)wyl, income (1 ? 而)wyl + ra, assets a and consumption c: 
% For the distributions, you may want to calculate Gini coefficients
% or if possible, plot the Lorenz curves, under the two different specifications. Is
% UBI welfare-improving? To answer this question you may want to compare the
% value functions v(a; y) under the two policies for some combinations of (a; y);
% or aggregate (utilitarian) social welfare

%% How I did it:
% I 


clear

%% Grid for labor - discrete labor choice
ifLabor = [0,1];
nGridLabor = length(ifLabor);

% households
rrho = 0.04; % discount rate
bbeta = 1/(1+rrho); % discount factor
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.95; % persistence of the income process
% ssigmaY = [0.2,0.8]; % income process variance
ssigmaY = 0.4; % income process variance
kkappa = 0.98465; % Calibrated in main_general_eq_getKappa.m

% firms
aalphaK = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
TFP = 1; % total production factor

% income shocks
nGridShocks = 31;
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
upperBound = 10 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?

% Grid of assets
a = 0.2;    % spacing of the grid (exponential)
nAssets = 100;

% Optimization
% options = optimset('Display','Iter','TolX',1e-07);  
options = optimset('Display', 'off','TolX',1e-04);

%% Government
nTao = 5;
vTao = linspace(0.125,0.25,nTao)';
llambda = 0.2;

vKOverL = zeros(nTao,1);
vDiff = zeros(nTao,1);
vGovBudgetGap = zeros(nTao,1);
vGovBudgetProfit = zeros(nTao,1);
vR = zeros(nTao,1);

for iTao = 1:nTao
    ttao = vTao(iTao);
    
    [kOverL,diff]=fminbnd(@(kOverL) ...
                        generalEqUBIFunction(kOverL,ttao,llambda, kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor),...
                        0.000001,10,options);

     [kOverLGap govBudgetGap govBudgetProfit r wage mValue mAssetPolicyIndex mAssetPolicy mConsumptionPolicy mLaborPolicy vStationaryDistribution capitalSupply laborSupplyEffective vIncomeShocks mTransition vGridAsset]=generalEqUBIFunction(kOverL,ttao,llambda, kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor);

    table(kOverL,diff,govBudgetGap,govBudgetProfit,r,wage)

    vKOverL(iTao) = kOverL;
    vDiff(iTao) = diff;
    vGovBudgetGap(iTao) = govBudgetGap;
    vGovBudgetProfit(iTao) = govBudgetProfit;
    
    vR(iTao) = r;
end


figure
plot(vTao,vKOverL)

figure
plot(vTao,vDiff)

figure
plot(vTao,vR)

figure
plot(vTao,vGovBudgetGap)

figure
plot(vTao,vGovBudgetProfit)

[v,ind] = min(vGovBudgetGap);

ttaoCalibrated = vTao(ind);

save('ttaoCalibrationManually')
                