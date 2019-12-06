%% Research Project 
% Economics 712
% Fall 2019
% Professor Dirk Krueger

% Value function iteration partially taken from "Basic RBC model with full depreciation" by Jesus Fernandez-Villaverde Haverford, July 3, 2013

% 2019-11-6 22:06:22
% 2019-12-3

clear;
% close all;

% Part II General Equilibrium - Computing Stationary Equilibria in the Aiyagari Model

% This part I write the Aiyagara model to MANUALLY calibrate how big the
% interest rate is, to get a sense of its degree.

%% model and parameters

% households
rrho = 0.04; % discount rate
bbeta = 1/(1+rrho); % discount factor
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.95; % persistence of the income process
% ssigmaY = [0.2,0.8]; % income process variance
ssigmaY = [0.2,0.4]; % income process variance

ssigmaError=ssigmaY*(sqrt(1-ddelta^2)); % variance of the error term of the income process
ssigmaErrorLow = min(ssigmaError);
ssigmaErrorHigh = max(ssigmaError);

% firms
aalphaK = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
labourEq = 1; % Normalize labor to 1
TFP = 1; % total production factor

% price
nR = 10;
% vGridR = linspace(-depreciation*(-0.2),rrho*1,nR)'; % interest rate range
% vGridR = linspace(0.02274,0.02277,nR)'; % interest rate range
vGridR = linspace(0.034239,0.034272,nR)'; % interest rate range

vCapitalGap = zeros(nR,1);
vCapitalSupply = zeros(nR,1);
vCapitalDemand = zeros(nR,1);

%% Grids for assets and shocks
% income shocks
nGridShocks = 31;
% **************************************************** TRY DIFFERENT VALUES ****************************************************
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorHigh,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
vIncomeShocks = exp(vIncomeShocks)';

% assets
AssetLimit = min(vIncomeShocks); % worst case scenario of income process
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
minK = 0;
% minK = -chi*AssetLimit/r; % natural borrowing limit

% **************************************************** TRY DIFFERENT VALUES ****************************************************
upperBound = 80 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?
maxK = upperBound; 
a = 0.2;    % spacing of the grid (exponential)
% **************************************************** TRY DIFFERENT VALUES ****************************************************
nAssets = 1000;
if chi == 0
    temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))';
    temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
    vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
%     clear temp temp_01;
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset));% life starting point of asset
    vGridAsset(iAsset0)=0;
else
    vGridAsset = curvspaceFunction(minK, maxK, nAssets ,1)';
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset)); % life starting point of asset
    vGridAsset(iAsset0)=0;
end

%% exercise 1 Infinite Horizon
mStationaryDistribution  = zeros(nGridAsset*nGridShocks,nR);

for iR = 1:nR
    r = vGridR(iR);    

    %% 1) Calculate demand for capital given r
    capitalDemand=(aalphaK*TFP/( r +depreciation))^(1/(1-aalphaK))*labourEq;

    %% 2) Calculate wage given r and demand for capital
    wage = (1-aalphaK)*TFP*(capitalDemand^aalphaK) * labourEq ^(-aalphaK); % wage

    %% 3) For a given interest rate guess, solve partial equilibrium to get the capital supply

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
                    consumption=(y* wage+(1+r)*asset-assetPrime);

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
                mConsumptionPolicy(iAsset,iShocks) = (y * wage + (1+r) * asset - assetPrime);
            end % asset

        end % shocks

        maxDifference = max(abs(mValueNew-mValue),[],'all');      
        mValue = mValueNew;
        iteration = iteration+1;

%         if mod(iteration,20)==1
%             fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration, maxDifference); 
%         end
    end
    toc
    fprintf(' Convergence achieved.\n Total Number of Iteration: %2.0f, Sup diff: %2.10f\n Interest rate: %2.6f\n', iteration, maxDifference, r); 
    
    % Obtain the total assets in the economy by calculating the
    % stationary distribution

    rows = zeros(1,nGridShocks*nGridAsset*nGridShocks);
    columns = zeros(1,nGridShocks*nGridAsset*nGridShocks);
    probs = zeros(1,nGridShocks*nGridAsset*nGridShocks);
    index = 1;

    for iAsset=1:nGridAsset % assets
        for iShocks=1:nGridShocks % shock today
            iAssetPrime = mAssetPolicyIndex(iAsset,iShocks);
            start = iAssetPrime;

            for iShocksPrime=1:nGridShocks % shock tomorrow
                rows(index) = nGridAsset*(iShocks-1)+iAsset;
                columns(index) = start+nGridAsset*(iShocksPrime-1);
                probs(index) = mTransition(iShocks,iShocksPrime);
                index=index+1;
            end
        end
    end

    ultimateTransition = sparse(rows,columns,probs,nGridShocks*nGridAsset,nGridShocks*nGridAsset);

    [vecP,~] = (eigs(ultimateTransition',1,'lr','Tolerance',1e-8));
    vecP = abs(vecP);
    vStationaryDistribution = vecP/sum(sum(vecP)); 
    mStationaryDistribution(:,iR)   = vStationaryDistribution;

    mPolFuncAssetVectorized=mAssetPolicy(:);
    capitalSupply=vStationaryDistribution'*mPolFuncAssetVectorized;
    
    vCapitalSupply(iR) = capitalSupply;
    vCapitalDemand(iR) = capitalDemand;   
    
    vCapitalGap(iR) = capitalSupply - capitalDemand ;
    fprintf('Supply - Demand = %2.6f\n',capitalSupply - capitalDemand)
end
figure 
plot(vGridR,vCapitalGap)

[v,iR_equilibrium] = min(abs(vCapitalGap));
R_equilibrium = vGridR(iR_equilibrium);

mStationaryDistributionEquilibrium=reshape(mStationaryDistribution(:,iR_equilibrium),[nGridAsset,nGridShocks]);
figure
[assetasset,shockshock]=meshgrid(log(vGridAsset), log(vIncomeShocks));
% [assetasset,shockshock]=meshgrid(vGridAsset, log(vIncomeShocks));
% linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets)
% [assetasset,shockshock]=meshgrid(exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))', log(vIncomeShocks));

mesh(assetasset, shockshock, mStationaryDistributionEquilibrium');
title('Stationary Equilibrium Distribution');
xlabel('Asset $log(a)$','interpreter','latex')
ylabel('Income Shocks $log(y)$','interpreter','latex')
zlabel('multiplier, in consumption units','interpreter','latex')
xlim([min(log(vGridAsset)),max(log(vGridAsset))])
ylim([min(log(vIncomeShocks)),max(log(vIncomeShocks))])
savefig('stationary_equilibrium_distribution')

save general_equilibrium

table(R_equilibrium)

% By plotting the distribution, we see the result is quite plausible. 
% since households are not allowed to borrow, a lot of them are stuck in
% the low income realm, but if we allow the household to borrow, the
% distribution would look much even.
% I use the large income shock variance, sigma = 0.4 to
% encourage the household to save. 

% Originally I got interest rate way bigger than rrho - the discount rate.
% Then I extended the grid of asset, and then things are normal finally.