clear
%% Question 3.3.2
% Evaluate Universal Basic Income

%% Step 2
% Then evaluate Andrew Young¡¯s proposal of a UBI: Compare the stationary
% equilibrium with ¦Ë = ¦Ó = 0 to a stationary equilibrium where ¦Ë = 0:2 and the
% tax rate ¦Ó adjusts to clear the government budget. 

%% How I did it:
% Given ttao, I used the fminbnd to find K/L
% Then in the outer loop, I used the fminbnd to find ttao to let government
% budget gap to be 0

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

% Government
nTao = 5;
vTao = linspace(0.125,0.25,nTao)';
llambda = 0.2;

tic
[ttao,gap]=fminbnd(@(ttao) ...
                    getTaoFunction(ttao,llambda, kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor,options),...
                    0.2188,0.25,options);
toc

[govBudgetGap,kOverL]=getTaoFunction(ttao,llambda, kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor,options);

 [kOverLGap govBudgetGap govBudgetProfit r wage mValue mAssetPolicyIndex mAssetPolicy mConsumptionPolicy mLaborPolicy vStationaryDistribution capitalSupply laborSupplyEffective vIncomeShocks mTransition vGridAsset] = generalEqUBIFunction(kOverL,ttao,llambda, kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor);

T = table(govBudgetGap,govBudgetProfit,r,wage,kkappa,ttao)

save('ttaoCalibration.mat')
save('withUBI')