% function [capitalGap laborParticipationRate] = generalEqFunction(r,ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,labourEq,TFP,...
%     nGridShocks,chi,upperBound,a,nAssets)
function [kOverLGap laborParticipationRate r wage mValue mAssetPolicyIndex mAssetPolicy mConsumptionPolicy mLaborPolicy vStationaryDistribution capitalSupply laborSupplyEffective vIncomeShocks mTransition vGridAsset] = generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor)
% function [kOverLGap] = generalEqUBIFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
%     nGridShocks,chi,upperBound,a,nAssets,...
%     ifLabor,nGridLabor)

ssigmaError=ssigmaY*(sqrt(1-ddelta^2)); % variance of the error term of the income process

%% Grids for assets and shocks
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaError,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
vIncomeShocks = exp(vIncomeShocks)';

AssetLimit = min(vIncomeShocks); % worst case scenario of income process
% minK = -chi*AssetLimit/r; % natural borrowing limit
minK = 0.0001;

maxK = upperBound; 

    vGridAsset = curvspaceFunction(minK, maxK, nAssets ,1)';
    nGridAsset = length(vGridAsset);


% if chi == 0
%     temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))';
%     temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
%     vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
% %     clear temp temp_01;
%     nGridAsset = length(vGridAsset);
%     [~, iAsset0]=min(abs(vGridAsset));% life starting point of asset
%     vGridAsset(iAsset0)=0;
% else
%     vGridAsset = curvspaceFunction(minK, maxK, nAssets ,1)';
%     nGridAsset = length(vGridAsset);
%     [~, iAsset0]=min(abs(vGridAsset)); % life starting point of asset
%     vGridAsset(iAsset0)=0;
% end

%% 1) Calculate r given kOverL
r = TFP * aalphaK * kOverL^(aalphaK-1) - depreciation;

%% 2) Calculate wage given kOverL
wage = (1-aalphaK)*TFP*(kOverL)^(aalphaK); % wage

%% 3) Given r and wage, solve partial equilibrium to get the supply side of kOverL 

%% Required Matrices and Vectors

% Initial Value Functions

mValue    = zeros(nGridAsset,nGridShocks); % using 0 as first Guess (average utility)
% mValueNew = zeros(nGridAsset,nGridShocks);
mValueNewTemp = zeros(nGridAsset,nGridShocks,nGridLabor);
mAssetPolicyTemp     = zeros(nGridAsset,nGridShocks,nGridLabor);
mConsumptionPolicyTemp = zeros(nGridAsset,nGridShocks,nGridLabor);
mAssetPolicyIndexTemp = zeros(nGridAsset,nGridShocks,nGridLabor);
% mLaborPolicy = zeros(nGridAsset,nGridShocks,nGridLabor);


%% Value function iteration

maxDifference = 10.0;
tolerance = 10^-8;
iteration = 0;
maxIter = 1000;

tic
while (maxDifference>tolerance) && iteration <= maxIter
    mExpectedValue = mValue*(mTransition');
    
    for iLabor = 1:nGridLabor
        labor = ifLabor(iLabor);

        for iShocks=1:nGridShocks
            y=vIncomeShocks(iShocks); 
            iAssetPrimeStart = 1; %including non-zero asset holding constraint

            for iAsset=1:nGridAsset
                asset=vGridAsset(iAsset);
                valueHighSoFar = -Inf;

%                 for iLabor = 1:nGridLabor
%                     labor = ifLabor(iLabor);                

                    for iAssetPrime = iAssetPrimeStart:nGridAsset
                        assetPrime=vGridAsset(iAssetPrime);
                        consumption=(y * wage * labor + (1+r)*asset - assetPrime);
                        
%                         if consumption <=0
%                             valueProvisional = - Inf;
                            
%                         else

                            valueProvisional = utilityFunction(consumption, ggamma) - labor * kkappa + ...
                                                                 bbeta * mExpectedValue(iAssetPrime,iShocks);

                            if (valueProvisional>valueHighSoFar) %&& consumption>0
                                valueHighSoFar = valueProvisional;
                                mAssetPolicyTemp(iAsset,iShocks,iLabor) = assetPrime;
                                iAssetPrimeStart = iAssetPrime;
                                mAssetPolicyIndexTemp(iAsset,iShocks,iLabor) = iAssetPrimeStart;
                            else
                                break; % We break when we have achieved the max
                            end
                            
%                         end% if consumption <0

                    end
                    mValueNewTemp(iAsset,iShocks,iLabor) = valueHighSoFar;
                    mConsumptionPolicyTemp(iAsset,iShocks,iLabor) = (y * wage * labor + (1+r) * asset - assetPrime);

%                 end%labor
  

            end % asset

        end % shocks
%     mLaborPolicy = mValueNewTemp(:,:,1)<mValueNewTemp(:,:,2);
    end%labor
    [mValueNew,mValueIndex] = max(mValueNewTemp,[],3); % max value on the third dimension
    mLaborPolicy = (mValueIndex~=1);
    
    mConsumptionPolicy = mConsumptionPolicyTemp(:,:,1) .* (mLaborPolicy == 0)  ...
        + mConsumptionPolicyTemp(:,:,2) .* (mLaborPolicy == 1);
    
    mAssetPolicy = mAssetPolicyTemp(:,:,1) .* (mLaborPolicy == 0)  ...
        + mAssetPolicyTemp(:,:,2) .* (mLaborPolicy == 1);
    
    mAssetPolicyIndex = mAssetPolicyIndexTemp(:,:,1) .* (mLaborPolicy == 0)  ...
        + mAssetPolicyIndexTemp(:,:,2) .* (mLaborPolicy == 1);
    
    maxDifference = max(abs(mValueNew-mValue),[],'all');      
    mValue = mValueNew;
    iteration = iteration+1;

        if mod(iteration,800)==1
            fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration, maxDifference); 
        end
end
toc
fprintf(' Convergence achieved.\n Total Number of Iteration: %2.0f, Sup diff: %2.10f\n Capital to Labor ratio: %2.6f\n', iteration, maxDifference, kOverL); 

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

mPolFuncAssetVectorized=mAssetPolicy(:);
capitalSupply=vStationaryDistribution'*mPolFuncAssetVectorized;

mLaborEffective = mLaborPolicy .* vIncomeShocks';
mLaborEffectiveVectorized = mLaborEffective(:);
laborSupplyEffective = vStationaryDistribution' * mLaborEffectiveVectorized;

kOverLSupply = capitalSupply/laborSupplyEffective;


kOverLGap = abs(kOverLSupply - kOverL) ;

laborParticipationRate = vStationaryDistribution' * (mLaborPolicy(:));

return