%% Research Project 
% Economics 712
% Fall 2019
% Professor Dirk Krueger

% Value function iteration partially taken from "Basic RBC model with full depreciation" by Jesus Fernandez-Villaverde Haverford, July 3, 2013

% 2019-11-6 22:06:22
% 2019-12-3

clear;
close all;

% Part I Partial Equilibrium

%% model and parameters
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.8; % persistence of the income process
ssigmaY = [0.2,0.4]; % income process variance
ssigmaError=ssigmaY*(sqrt(1-ddelta^2)); % variance of the error term of the income process
ssigmaErrorLow = min(ssigmaError);
ssigmaErrorHigh = max(ssigmaError);

rrho = 0.04; % discount rate
bbeta = 1/(1+rrho); % discount factor
r = 0.02;% interest rate

% income shocks
nGridShocks = 11;
% **************************************************** TRY DIFFERENT VALUES ****************************************************
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorLow,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
vIncomeShocks = exp(vIncomeShocks)';

% assets
AssetLimit = min(vIncomeShocks); % worst case scenario of income process
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
minK = -chi*AssetLimit/r; % natural borrowing limit

% **************************************************** TRY DIFFERENT VALUES ****************************************************
upperBound = 4 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?
maxK = upperBound; 
a = 0.2;    % spacing of the grid (exponential)
% **************************************************** TRY DIFFERENT VALUES ****************************************************
nAssets = 500;
temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))';
temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
clear temp temp_01;
nGridAsset = length(vGridAsset);
[~, index]=min(abs(vGridAsset));
vGridAsset(index)=0;

%% exercise 1
% Formulate the problem of the agent recursively, i.e. write down Bellman¡¯s equation 
% and derive the stochastic Euler equation.
% Bellman¡¯s equation
% v(a,y) = max @a' { (1-bbeta)*u(c)+bbeta*Ev(a',y')} s.t. c+a' = y+(1+r)*a, and a'>=0
% stochastic Euler equation
% c^(-ggamma) = bbeta * (1 + r) * E(c' ^(-ggama)) + mmu/(1-bbeta); % Note mmu>=0 is the lagrangian multiplier of the borrowing constraint a'>=o

%% exercise 2 & 4 - infinite horizon and simulation
nSimulations = 2000; % number of simulated paths M
nPeriods = 61; % number of periods in one simulation

%% Required Matrices and Vectors

% Initial Value Functions

mValue    = zeros(nGridAsset,nGridShocks); % using 0 as first Guess (average utility)
mValueNew = zeros(nGridAsset,nGridShocks);
mAssetPolicy     = zeros(nGridAsset,nGridShocks);
mConsumptionPolicy = zeros(nGridAsset,nGridShocks);
mAssetPolicyIndex = zeros(nGridAsset,nGridShocks);

%% Value function iteration

maxDifference = 10.0;
tolerance = 10^-8;
iteration = 0;

tic
while (maxDifference>tolerance)
    mExpectedValue = mValue*(mTransition');

    for iShocks=1:nGridShocks
        y=vIncomeShocks(iShocks); 
        iAssetPrimeStart = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset
            asset=vGridAsset(iAsset);
            valueHighSoFar = -Inf;

            for iAssetPrime = iAssetPrimeStart:nGridAsset
                assetPrime=vGridAsset(iAssetPrime);
                consumption=(y+(1+r)*asset-assetPrime);

                valueProvisional = (1-bbeta) * utilityFunction(consumption, ggamma) + ...
                                                     bbeta * mExpectedValue(iAssetPrime,iShocks);

                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    mAssetPolicy(iAsset,iShocks) = assetPrime;
                    iAssetPrimeStart = iAssetPrime;
                    mAssetPolicyIndex(iAsset,iShocks) = iAssetPrimeStart;
                else
                    break; % We break when we have achieved the max
                end    
            end
                mValueNew(iAsset,iShocks) = valueHighSoFar;
                mConsumptionPolicy(iAsset,iShocks) = (y + (1+r) * asset - assetPrime);
        end % asset

    end % shocks
    
    maxDifference = max(abs(mValueNew-mValue),[],'all');      
    mValue = mValueNew;
    iteration = iteration+1;
    
    if mod(iteration,20)==1
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration, maxDifference); 
    end
end
toc
fprintf(' Convergence achieved.\n Total Number of Iteration: %2.0f, Sup diff: %2.10f\n', iteration, maxDifference); 

% lagrangian multiplier of borrowing constraint in consumption units
lagMultiplierBorrowingConstraintInConsumptionUnits = lagMultiplierBorrowingConstraintFunction(mAssetPolicyIndex,mConsumptionPolicy, ggamma, nGridAsset, nGridShocks, vIncomeShocks, mTransition, vGridAsset, bbeta, r);
% some of the multipliers are negative. WHY???

%% simulation
mShocksIndexSimulation= ash_panelFunction(mTransition,nSimulations,nPeriods,nGridShocks)'; % nPeriods by nSimulations
mAssetIndexSimulation = ones(nPeriods,nSimulations); % Note this is an index matrix, so start with 1
mConsumptionSimulation = zeros(nPeriods,nSimulations);

for iSimulations = 1 : nSimulations
    iAsset = 1;
    for iPeriods = 1 : nPeriods
        iShocks = mShocksIndexSimulation(iPeriods,iSimulations);
        mAssetIndexSimulation(iPeriods+1,iSimulations) = mAssetPolicyIndex(iAsset,iShocks);
        mConsumptionSimulation(iPeriods,iSimulations) = mConsumptionPolicy(iAsset,iShocks);
        iAsset = mAssetPolicyIndex(iAsset,iShocks);
    end
end

mAssetIndexSimulation = mAssetIndexSimulation(1:nPeriods,:);

save('infinite_horizon_low_shock')


%% figures

[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);

figure

subplot(1,2,1)
load infinite_horizon_low_shock;
mesh(assetasset, shockshock, lagMultiplierBorrowingConstraintInConsumptionUnits');
title('Lag Multiplier of Borrowing Constraint - Infinite Horizon - low variance shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('multiplier, in consumption units','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load infinite_horizon_high_shock;
mesh(assetasset, shockshock, lagMultiplierBorrowingConstraintInConsumptionUnits');
title('Lag Multiplier of Borrowing Constraint - Infinite Horizon - high variance shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('lag_multiplier_infinite_horizon')


figure
subplot(1,2,1)
load infinite_horizon_low_shock;
mesh(assetasset, shockshock, mValue');
title('Value - Infinite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('value','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load infinite_horizon_high_shock;
mesh(assetasset, shockshock, mValue');
title('Value - Infinite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('value','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('value_infinite_horizon')


figure
subplot(1,2,1)
load infinite_horizon_low_shock;
mesh(assetasset, shockshock, mAssetPolicy');
title('Policy for Next Period Asset Holdings - Infinite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load infinite_horizon_high_shock;
mesh(assetasset, shockshock, mAssetPolicy');
title('Policy for Next Period Asset Holdings - Infinite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('assetPolicy_infinite_horizon')


figure
subplot(1,2,1)
load infinite_horizon_low_shock;
mesh(assetasset, shockshock, mConsumptionPolicy');
title('Policy for Consumption - Infinite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('consumption','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load infinite_horizon_high_shock;
mesh(assetasset, shockshock, mConsumptionPolicy');
title('Policy for Consumption - Infinite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('consumption','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('consumptionPolicy_infinite_horizon')


figure 
load infinite_horizon_low_shock;

subplot(2,2,1)
plot(1:nPeriods,mean(mConsumptionSimulation,2))
xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Consumption - Low Variance Shocks')
xlabel('T')
ylabel('c')

subplot(2,2,2)
plot(1:nPeriods,mean(vGridAsset(mAssetIndexSimulation),2))
xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Asset Holdings - Low Variance Shocks')
xlabel('T')
ylabel('a')

load infinite_horizon_high_shock;
subplot(2,2,3)
plot(1:nPeriods,mean(mConsumptionSimulation,2))
xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Consumption - High Variance Shocks')
xlabel('T')
ylabel('c')

subplot(2,2,4)
plot(1:nPeriods,mean(vGridAsset(mAssetIndexSimulation),2))
xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Asset Holdings - High Variance Shocks')
xlabel('T')
ylabel('a')

savefig('simulation_infinite_horizon')

%% plot consumption function and interpret the difference

figure
load('infinite_horizon_low_shock.mat')
plot(vGridAsset,mConsumptionPolicy,'b');

hold on
load('infinite_horizon_high_shock.mat')
plot(vGridAsset,mConsumptionPolicy,'r');

xlabel('asset');
ylabel('consumption');
title('Consumption Under Low (blue) and High (red) Shocks')
savefig('consumptionPolicy_comparison_infinite_horizon_2D')

% We can see from the consumption plot that 
% under high income shocks, households consume more and vice versa.