%% Research Project 
% Economics 712
% Fall 2019
% Professor Dirk Krueger

% Value function iteration partially taken from "Basic RBC model with full depreciation" by Jesus Fernandez-Villaverde Haverford, July 3, 2013

% 2019-11-6 22:06:22
% 2019-12-3

clear;
% close all;

%% Added income profile, income replacement, and mortality
theta=0.7;  %income replacement ratio
load incprofile.txt   %45-by-1
load survs.txt        %61-by-1 
load consprofile.txt  %264-by-2


% Part I Partial Equilibrium

%% model and parameters
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.95; % persistence of the income process
ssigmaY = [0.2,0.25]; % income process variance
ssigmaError=ssigmaY*(sqrt(1-ddelta^2)); % variance of the error term of the income process
ssigmaErrorLow = min(ssigmaError);
ssigmaErrorHigh = max(ssigmaError);

rrho = 0.03; % discount rate
bbeta = 1/(1+rrho); % discount factor
r = 0.02;% interest rate

% income shocks
nGridShocks = 11;
% **************************************************** TRY DIFFERENT VALUES ****************************************************
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorHigh,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
vIncomeShocks = exp(vIncomeShocks)';

% assets
AssetLimit = min(vIncomeShocks); % worst case scenario of income process
chi = 0.05; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
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
    vGridAsset = curvspaceFunction(minK, maxK, nAssets ,2)';
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset)); % life starting point of asset
    vGridAsset(iAsset0)=0;
end

%% exercise 7 & 8 & 9 Finite Horizon with retirement and death risk
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

tic
for iPeriods = nPeriods:-1:1
    psurv = survs(iPeriods);

    mExpectedValue = mValue(:,:,iPeriods+1)*(mTransition');
    for iShocks=1:nGridShocks
        if iPeriods>45
            y=theta*incprofile(end)*vIncomeShocks(iShocks); 
        else
            y=incprofile(iPeriods)*vIncomeShocks(iShocks); 
        end
%                     y=vIncomeShocks(iShocks); 

        
        iAssetPrimeStart = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset
            asset=vGridAsset(iAsset);
            valueHighSoFar = -Inf;

            for iAssetPrime = iAssetPrimeStart:nGridAsset
                assetPrime=vGridAsset(iAssetPrime);
                consumption=(y+(1+r)*asset-assetPrime);
                valueProvisional = utilityFunction(consumption, ggamma) + ...
                                     psurv*bbeta * mExpectedValue(iAssetPrime,iShocks);

                if (valueProvisional>valueHighSoFar)
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

age = 19+(1:nPeriods);
agedata  = consprofile(:,1);
consdata = consprofile(:,2);
figure
plot(age,mean(mConsumptionSimulation'),agedata,consdata,'--','linewidth',1.5)
legend('model','data','location','northeast')

% The levels of the two profiles are different become income is scaled
% differently, so first of all we need to rescale incprofile to get the
% level right. Next, changing the discount rate and the borrowing limit
% can improve the fit a lot.

%% Q9 consumption insurance coeffcient, by age
phis = zeros(nPeriods-1,1);
for tt=2:nPeriods
    dlogc = log(mConsumptionSimulation(tt,:))'-log(mConsumptionSimulation(tt-1,:))';
    dlogy = log(vIncomeShocks(mShocksIndexSimulation(tt,:)))-log(vIncomeShocks(mShocksIndexSimulation(tt-1,:)));
    covcy = cov(dlogc,dlogy);  %this is 2-by-2 variance-covariance matrix
    phis(tt-1) = 1-covcy(1,2)/covcy(2,2);
end
age = 19+(2:nPeriods);
% figure
% plot(age,phis,'linewidth',1.5)
% title(['consumption insurance coeffcient, \delta=',num2str(ddelta)])

% save( 'ddeltaEqual0','phis')
% save( 'ddeltaEqual095','phis')
% save( 'ddeltaEqual099','phis')

figure
load ddeltaEqual0
plot(age,phis,'r','linewidth',1.5)
hold on
load ddeltaEqual095
plot(age,phis,'b','linewidth',1.5)
hold on
load ddeltaEqual099
plot(age,phis,'g','linewidth',1.5)
title(['consumption insurance coeffcient, \delta=',num2str(ddelta)])
legend('$\delta = 0$','$\delta = 0.95$','$\delta = 0.99$','data','location','northeast','interpreter','latex')
savefig('consumption_insurance_under_different_persistence')

% By changing ddelta, the persistence coefficient, we can see that as
% persistence grows, the coefficient goes down. When there is no
% persistence at all, ddelta = 0, the coefficient is almost 1, and when
% there is perfect persistence, ddelta = 1, or almost, the coefficient is
% almost 0. Such is intuitive, since if income process is very much
% persistent, consumption would covary to a very large extent with income,
% but when there is iid income process, consumption risk can be perfectly
% insured away, so the coefficient is almost 1.