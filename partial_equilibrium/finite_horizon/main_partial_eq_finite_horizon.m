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
ddelta = 0.95; % persistence of the income process
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
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorHigh,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
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
if chi == 0
    temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))';
    temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
    vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
    clear temp temp_01;
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset));% life starting point of asset
    vGridAsset(iAsset0)=0;
else
    vGridAsset = curvspaceFunction(minK, maxK, nAssets ,1)';
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset)); % life starting point of asset
    vGridAsset(iAsset0)=0;
end

%% exercise 3 & 5 & 6 Finite Horizon
nSimulations = 50000; % number of simulated paths M
nPeriods = 61; % horizon

%% Required Matrices and Vectors

% Initial Value Functions

mValue    = zeros(nGridAsset,nGridShocks,nPeriods+1); % using 0 as final period (T+1) value
if chi ~= 0
    mValueDeath = zeros(nGridAsset,nGridShocks);
    mValueDeath(1:iAsset0-1,:) = -100000000;%*ones(iAsset0-1,nGridShocks);
    mValue(:,:,end) = mValueDeath;
end
mAssetPolicy     = zeros(nGridAsset,nGridShocks,nPeriods);
mConsumptionPolicy = zeros(nGridAsset,nGridShocks,nPeriods);
mAssetPolicyIndex = ones(nGridAsset,nGridShocks,nPeriods); % index matrix starts from all 1s

%% Value function iteration
% We can do grid search and vectorize eveything by not iterating over assetPrime 
% but using vGridAsset as assetPrime, 
% but apart from being considerably slower than the following method, we
% also have to take care of the complex numbers since our utility function
% involves exponentials.
% 
% tic
% for iPeriods = nPeriods:-1:1
%     mExpectedValue = mValue(:,:,iPeriods+1)*(mTransition');
% 
%     for iShocks=1:nGridShocks
%         y=vIncomeShocks(iShocks); 
%         iAssetPrimeStart = 1; %including non-zero asset holding constraint
% 
%         for iAsset=1:nGridAsset
%             asset=vGridAsset(iAsset);
%             valueHighSoFar = -Inf;
% 
%             for iAssetPrime = iAssetPrimeStart:nGridAsset
%                 assetPrime=vGridAsset(iAssetPrime);
%                 consumption=(y+(1+r)*asset-assetPrime);
%                 valueProvisional = utilityFunction(consumption, ggamma) + ...
%                                      bbeta * mExpectedValue(iAssetPrime,iShocks);
% %                  if valueProvisional == -Inf && valueHighSoFar == -Inf
% %                  if (valueProvisional == -Inf && valueHighSoFar == -Inf) | consumption <= 0
%                 if (valueProvisional>valueHighSoFar)
% %                     valueHighSoFar = valueProvisional;
%                      
% %                  elseif (valueProvisional>valueHighSoFar)
%                     valueHighSoFar = valueProvisional;
%                     mAssetPolicy(iAsset,iShocks,iPeriods) = assetPrime;
%                     iAssetPrimeStart = iAssetPrime;
%                     mAssetPolicyIndex(iAsset,iShocks,iPeriods) = iAssetPrimeStart;                    
%                  else
%                     break; % We break when we have achieved the max
%                 end                       
%             end
%             mValue(iAsset,iShocks,iPeriods) = valueHighSoFar;
%             mConsumptionPolicy(iAsset,iShocks,iPeriods) = (y + (1+r) * asset - assetPrime);            
% 
%         end
%     end
% 
% %     if mod(iPeriods,20)==0
%         fprintf(' Time Periods: %2.0f \n', iPeriods); 
% %     end
% end
% toc
% fprintf('Iteration finished.\n'); 
% 
% mValue = mValue(:,:,1:nPeriods);


tic
for iPeriods = nPeriods:-1:1
    if iPeriods == nPeriods
        iAssetPrime =iAsset0;
        mAssetPolicyIndex(:,:,iPeriods) = iAsset0;  
        mAssetPolicy(:,:,iPeriods) = 0;  
        mConsumptionPolicy(:,:,iPeriods) = vIncomeShocks' + (1+r)* vGridAsset;
        for iShocks=1:nGridShocks
            for iAsset=1:nGridAsset
                consumption = mConsumptionPolicy(iAsset,iShocks,iPeriods);
                     
                mValue(iAsset,iShocks,iPeriods) = utilityFunction(consumption, ggamma);
            end
        end
                    
    else

        mExpectedValue = mValue(:,:,iPeriods+1)*(mTransition');

        for iShocks=1:nGridShocks
            y=vIncomeShocks(iShocks); 
            iAssetPrimeStart = 1; %including non-zero asset holding constraint

            for iAsset=1:nGridAsset
                asset=vGridAsset(iAsset);
                valueHighSoFar = -Inf;

                for iAssetPrime = iAssetPrimeStart:nGridAsset
                    assetPrime=vGridAsset(iAssetPrime);
                    consumption=(y+(1+r)*asset-assetPrime);
                    valueProvisional = utilityFunction(consumption, ggamma) + ...
                                         bbeta * mExpectedValue(iAssetPrime,iShocks);
    %                  if valueProvisional == -Inf && valueHighSoFar == -Inf
    %                  if (valueProvisional == -Inf && valueHighSoFar == -Inf) | consumption <= 0
                    if (valueProvisional>valueHighSoFar)
    %                     valueHighSoFar = valueProvisional;

    %                  elseif (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mAssetPolicy(iAsset,iShocks,iPeriods) = assetPrime;
                        iAssetPrimeStart = iAssetPrime;
                        mAssetPolicyIndex(iAsset,iShocks,iPeriods) = iAssetPrimeStart;                    
                     else
                        break; % We break when we have achieved the max
                    end                       
                end
                mValue(iAsset,iShocks,iPeriods) = valueHighSoFar;
                mConsumptionPolicy(iAsset,iShocks,iPeriods) = (y + (1+r) * asset - assetPrime);            

            end
        end

    %     if mod(iPeriods,20)==0
        fprintf(' Time Periods: %2.0f \n', iPeriods); 
    %     end
    end
end
toc
fprintf('Iteration finished.\n'); 

mValue = mValue(:,:,1:nPeriods);

%% simulation
mShocksIndexSimulation= ash_panelFunction(mTransition,nSimulations,nPeriods,nGridShocks)'; % nPeriods by nSimulations
mAssetHoldingsIndexSimulation =[ iAsset0 * ones(1,nSimulations);ones(nPeriods,nSimulations)]; % Note this is an index matrix, so start with 1
mConsumptionSimulation = zeros(nPeriods,nSimulations);

for iSimulations = 1 : nSimulations
    iAsset = iAsset0;
    for iPeriods = 1 : nPeriods
        iShocks = mShocksIndexSimulation(iPeriods,iSimulations);
        mAssetHoldingsIndexSimulation(iPeriods+1,iSimulations) = mAssetPolicyIndex(iAsset,iShocks,iPeriods);                     
        mConsumptionSimulation(iPeriods,iSimulations) = mConsumptionPolicy(iAsset,iShocks,iPeriods);
        iAsset = mAssetPolicyIndex(iAsset,iShocks,iPeriods);
    end
end

mAssetHoldingsIndexSimulation = mAssetHoldingsIndexSimulation(1:nPeriods,:);
save('finite_horizon_high_shock')


%% figures

figure
subplot(1,2,1)
load finite_horizon_low_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mValue(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mValue(:,:,iPeriods)');
end

title('Value - Finite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('value','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load finite_horizon_high_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mValue(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mValue(:,:,iPeriods)');
end
title('Value - Finite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('value','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])
% zlim([-5,5])
savefig('value_finite_horizon')


figure
subplot(1,2,1)
load finite_horizon_low_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mAssetPolicy(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mAssetPolicy(:,:,iPeriods)');
end
title('Policy for Next Period Asset Holdings - Finite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load finite_horizon_high_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mAssetPolicy(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mAssetPolicy(:,:,iPeriods)');
end
title('Policy for Next Period Asset Holdings - Finite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('assetPolicy_finite_horizon')


figure
subplot(1,2,1)
load finite_horizon_low_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mConsumptionPolicy(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mConsumptionPolicy(:,:,iPeriods)');
end

title('Policy for Consumption - Finite Horizon - Low Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('consumption','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

subplot(1,2,2)
load finite_horizon_high_shock;
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
mesh(assetasset, shockshock, mConsumptionPolicy(:,:,1)');

for iPeriods = 2:nPeriods
    hold on
    mesh(assetasset, shockshock, mConsumptionPolicy(:,:,iPeriods)');
end
title('Policy for Consumption - Finite Horizon - High Variance Shocks','interpreter','latex')
xlabel('Asset $a$','interpreter','latex')
ylabel('Income Shocks $y$','interpreter','latex')
zlabel('consumption','interpreter','latex')
xlim([min(vGridAsset),max(vGridAsset)])
ylim([min(vIncomeShocks),max(vIncomeShocks)])

savefig('consumptionPolicy_finite_horizon')


figure 
subplot(1,2,1)
load finite_horizon_low_shock;
plot(1:nPeriods,mean(mConsumptionSimulation,2),'b')
hold on 
load finite_horizon_high_shock;
plot(1:nPeriods,mean(mConsumptionSimulation,2),'r')
xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Consumption under Different Shocks')
legend('low variance shocks','high variance shocks')
xlabel('T')
ylabel('c')

subplot(1,2,2)
load finite_horizon_low_shock;
plot(1:nPeriods,mean(vGridAsset(mAssetHoldingsIndexSimulation),2),'b')
hold on
load finite_horizon_high_shock;
plot(1:nPeriods,mean(vGridAsset(mAssetHoldingsIndexSimulation),2),'r')

xlim([1,nPeriods]);
% ylim([0,max(max(mean(mConsumptionSimulation,2)),max(mean(vGridAsset(mAssetIndexSimulation),2)))]);
title('Average Asset Holdings under Different Shocks')
legend('low variance shocks','high variance shocks')
xlabel('T')
ylabel('a')

savefig('simulation_finite_horizon')

%% plot consumption function and interpret the difference

figure
load('finite_horizon_low_shock.mat')
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
plot(vGridAsset,mConsumptionPolicy(:,:,1),'b');

hold on
load('finite_horizon_high_shock.mat')
[assetasset,shockshock]=meshgrid(vGridAsset, vIncomeShocks);
plot(vGridAsset,mConsumptionPolicy(:,:,56),'r');

xlabel('asset');
ylabel('consumption');
title('Consumption of Young (blue T=0) and Old (red T=55) Households')
savefig('consumptionPolicy_comparison_finite_horizon_2D')

% We can see from the graph that unless the shocks are really really bad,
% the consumption of older households are bigger than younger ones.

open simulation_finite_horizon.fig
% I didn't get a hump shape for the lifetime consumption profile. The
% consumption profile seems to increase over life time. 
% Maybe it's because the household is too patient to spend during the prime
% of its life. So rrho should be larger, or interest rate r be lower, or
% households are allowed to borrow.

% I tried everything but still can't get the hump shape.
% Maybe it's because the household always has income risks, unlike
% retirement benefit.

% Once I added the borrowing constraint, with the original persistance
% parameter, the asset holding files are always U-shaped, exactly the
% opposite of hump shape. And the consumption goes down persistently 
% throughout life. So I make the household more patient by lowering
% rrho below the interest rate so that they want to save. I also decrease
% the income shock variance. But none of this works.

%================Reason for Faiilure===============
% I can't get the hump shape.
% I now get the two inferences up till now:
% If the households are ALLOWED to borrow, the problem lies in the fact that
% the interst rate is exogenous. More spesicifically, the value function is negative
% infinity for the most entry of the value faction, 3D, and the reason is
% that the interest rate is exogenous, so for some states of the last
% periods, the value would always be negative infinity, or very very negative,
% and that would cause the whole value function behaves wierd.

% If the households are NOT allowed to borrow, the problem lies in the fact
% that the households will always face income shock, even in the last
% period of his life. That would cause the consumption to accumulate till
% the very last periods, instead of U-shaped like in the empirical data,
% for which the income after retirement is quite deterministic compared to
% one's salaries. Such can be seen once we did question 7-9, where we can
% get the consumption profile to be hump-shaped.