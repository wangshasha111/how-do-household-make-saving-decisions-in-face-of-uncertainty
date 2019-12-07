%% Question 3.3.1
% Evaluate Universal Basic Income
% I search manually for kkappa to calibrate labor paticipation to be 0.8
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

% kkappa = 0.5;
nKappa = 10;
vKappa = linspace(0.9111,1,nKappa)';
vKOverL = zeros(nKappa,1);
vDiff = zeros(nKappa,1);
vLaborParticipationRate = zeros(nKappa,1);

%% Given kkappa, solve for the general equilibrium

for iKappa = 1:nKappa
    kkappa = vKappa(iKappa);

    [kOverL,diff]=fminbnd(@(kOverL) ...
                        generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
                        nGridShocks,chi,upperBound,a,nAssets,...
                        ifLabor,nGridLabor),...
                        0.000001,10,options);

    [kOverLGap laborParticipationRate r wage mValue mAssetPolicyIndex mAssetPolicy mConsumptionPolicy mLaborPolicy vStationaryDistribution capitalSupply laborSupplyEffective  vIncomeShocks mTransition vGridAsset] =generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
        nGridShocks,chi,upperBound,a,nAssets,...
        ifLabor,nGridLabor)

    table(kOverL,diff,laborParticipationRate)

    vKOverL(iKappa) = kOverL;
    vDiff(iKappa) = diff;
    vLaborParticipationRate(iKappa) = laborParticipationRate;
                
end

figure
plot(vKappa,vKOverL)

figure
plot(vKappa,vDiff)

figure
plot(vKappa,vLaborParticipationRate)

[v,ind] = min(vLaborParticipationRate-0.8);

kkappaCalibrated = vKappa(ind);

save('kkappaCalibration')