clear

%% Question 3.3.3 Policy Evaluation
% What happens to macroeconomic aggregates (Y; K; C), to equilibrium prices (w; r) 
% and the equilibrium distributions for earnings (1 -¦Ó)wyl, income (1 - ¦Ó)wyl + ra, assets a
% and consumption c.

% For the distributions, you may want to calculate Gini coefficients
% or if possible, plot the Lorenz curves, under the two different specifications. 

% Is UBI welfare-improving? To answer this question you may want to compare the
% value functions v(a; y) under the two policies for some combinations of (a; y);
% or aggregate (utilitarian) social welfare

load('withoutUBI.mat')

vWage = zeros(1,2);
vRate = zeros(1,2);
vAggY = zeros(1,2);
vAggK = zeros(1,2);
vAggConsumption = zeros(1,2);

vDistribution = zeros(nGridShocks*nAssets,2);
vIndEarnings = zeros(nGridShocks*nAssets,2);
vIndIncome = zeros(nGridShocks*nAssets,2);
vIndAsset = zeros(nGridShocks*nAssets,2);
vIndConsumption = zeros(nGridShocks*nAssets,2);
vGiniCoefficients = zeros(1,2);

vDistribution(:,1) = vStationaryDistribution;

vWage(1) = wage;
vRate(1) = r;
vAggK(1) = capitalSupply;
vAggY(1) = capitalSupply^aalphaK * laborSupplyEffective^(1-aalphaK);

vConsumption = mConsumptionPolicy(:);
vAggConsumption(1) = vStationaryDistribution' * vConsumption;

mIndEarnings = mLaborPolicy.*vInc
vIndEarnings(:,1) = 
vIndIncome(:,1) = 
vIndAsset(:,1) = 
vIndConsumption(:,1) = 
vGiniCoefficients(1) = 

%==========================================
load('withUBI.mat')
vDistribution(:,2) = vStationaryDistribution;

vWage(2) = wage;
vRate(2) = r;
vAggK(2) = capitalSupply;
vAggY(2) = capitalSupply^aalphaK * laborSupplyEffective^(1-aalphaK);

vConsumption = mConsumptionPolicy(:);
vAggConsumption(2) = vStationaryDistribution' * vConsumption;

vIndEarnings(:,2) = 
vIndIncome(:,2) = 
vIndAsset(:,2) = 
vIndConsumption(:,2) = 
vGiniCoefficients(2) = 