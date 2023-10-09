%% 1.  Transition without UBI % I failed to get it converge
% Again choose T = 1 and your preferred parameterization from the last
% question. Compute the transition path induced by an unexpected Great
% Recession. Specifically, in the economy without UBI, start with a stationary 
% equilibrium in which aggregate TFP is equal to 1. Now suppose that
% for 10 years TFP falls by 10 % and then recovers back to 1. This event
% is perfectly unexpected, but once it occurs the path of TFP is perfectly
% foreseen. That is, in the production function

% Yt = At F(Kt; Lt)

% the productivity term equals At = 1 for t = 0 and t > 11; and At = 0:9
% for 1 <= t <= 10: The change in TFP induces a transition path from the
% initial steady state to the final steady state (which is equal to the initial
% steady state) that should look like a recession. Note that you have already
% computed the initial (final) steady state, so you know where the transition
% starts and where it ends.

% Plot the time paths of aggregate TFP, output, capital and consumption 
% (in percentage deviations from the initial steady state). 
% How could you alter the model to generate more amplification of
% the TFP shocks?

clear;

load withoutUBI.mat

% 2) Guess a sequence of aggregate capital
kOverLInitial = kOverL;
kOverLFinal = kOverL;
nPeriods = 60; % T
vKOverLPathGuess = linspace(kOverLInitial,kOverLFinal,nPeriods)';

vDistributionInitial         = vStationaryDistribution;
vDistributionFinal           = vStationaryDistribution;
mValueInitial       = mValue;
mValueFinal = mValue;
mConsumptionPolicyInitial = mConsumptionPolicy;
mConsumptionPolicyFinal   = mConsumptionPolicy;
capitalSupplyInitial       = capitalSupply;
capitalSupplyFinal         = capitalSupply;
mAssetPolicyInitial       = mAssetPolicy;
mAssetPolicyFinal         = mAssetPolicy;
mAssetPolicyIndexInitial       = mAssetPolicyIndex;
mAssetPolicyIndexFinal         = mAssetPolicyIndex;
mLaborPolicyInitial = mLaborPolicy;
mLaborPolicyFinal = mLaborPolicy;


% 3) Given initial guess, calculate interest rates and
% find new value functions backwards. Given initial distribution,
% Iterate forward to find new guess. Use relaxation.

mValuePath    = zeros(nAssets,nGridShocks,nPeriods); 
mValuePath(:,:,end) = mValueFinal;

mAssetPolicyPath     = zeros(nAssets,nGridShocks,nPeriods);
mAssetPolicyPath(:,:,end) = mAssetPolicyFinal;

mAssetPolicyIndexPath = ones(nAssets,nGridShocks,nPeriods); % index matrix starts from all 1s
mAssetPolicyIndexPath(:,:,end) = mAssetPolicyIndexFinal;

mConsumptionPolicyPath = zeros(nAssets,nGridShocks,nPeriods);
mConsumptionPolicyPath(:,:,end) = mConsumptionPolicyFinal;

mLaborPolicyPath = zeros(nAssets,nGridShocks,nPeriods);
mLaborPolicyPath(:,:,end) = mLaborPolicyFinal;

mValueNewTemp = zeros(nAssets,nGridShocks,nGridLabor);
mAssetPolicyTemp     = zeros(nAssets,nGridShocks,nGridLabor);
mConsumptionPolicyTemp = zeros(nAssets,nGridShocks,nGridLabor);
mAssetPolicyIndexTemp = zeros(nAssets,nGridShocks,nGridLabor); % These are the matrices for labor choice comparison

vKOverLPathNew = vKOverLPathGuess;

%% Find transition path

% distanceOffStability = 100;
tol = 1e-03;
iteration = 1;
iterMax = 100;
vDistanceOffStability = zeros(iterMax+1,1);
vDistanceOffStability(1) = 100;

while vDistanceOffStability(iteration) >= tol && iteration <= iterMax

    vTFP = ones(nPeriods,1);    
    vTFP(1:10)=0.9;
    vR = vTFP * aalphaK .* vKOverLPathNew.^(aalphaK-1) - depreciation;
    vWage = vTFP * (1 - aalphaK) .* vKOverLPathNew.^(aalphaK);

    % Value function iteration

    tic
    for iPeriods = nPeriods-1:-1:1

        wage = vWage(iPeriods);
        r = vR(iPeriods);
        mExpectedValue = mValuePath(:,:,iPeriods+1)*(mTransition');

        for iLabor = 1:nGridLabor
            labor = ifLabor(iLabor);


            for iShocks=1:nGridShocks
                y=vIncomeShocks(iShocks); 
                iAssetPrimeStart = 1; %including non-zero asset holding constraint

                for iAsset=1:nAssets
                    asset=vGridAsset(iAsset);

                    valueHighSoFar = -Inf;

                    for iAssetPrime = iAssetPrimeStart:nAssets
                        assetPrime=vGridAsset(iAssetPrime);
                        consumption=(y * wage * labor + (1 + r) * asset - assetPrime);

                        valueProvisional = utilityFunction(consumption, ggamma)  - labor * kkappa + ...
                                             bbeta * mExpectedValue(iAssetPrime,iShocks);

                        if (valueProvisional>valueHighSoFar)
                            valueHighSoFar = valueProvisional;
                            mAssetPolicyTemp(iAsset,iShocks,iLabor) = assetPrime;
                            iAssetPrimeStart = iAssetPrime;
                            mAssetPolicyIndexTemp(iAsset,iShocks,iLabor) = iAssetPrimeStart;
                         else
                            break; % We break when we have achieved the max
                        end                       
                    end % assetPrime
                    mValueNewTemp(iAsset,iShocks,iLabor) = valueHighSoFar;
                    mConsumptionPolicyTemp(iAsset,iShocks,iLabor) = (y * wage * labor + (1+r) * asset - assetPrime);

                end%asset
            end%shocks
        end%labor

        [mValueNew,mValueIndex] = max(mValueNewTemp,[],3); % max value on the third dimension
        mLaborPolicy = (mValueIndex~=1);

        mConsumptionPolicy = mConsumptionPolicyTemp(:,:,1) .* (mLaborPolicy == 0)  ...
            + mConsumptionPolicyTemp(:,:,2) .* (mLaborPolicy == 1);

        mAssetPolicy = mAssetPolicyTemp(:,:,1) .* (mLaborPolicy == 0)  ...
            + mAssetPolicyTemp(:,:,2) .* (mLaborPolicy == 1);

        mAssetPolicyIndex = mAssetPolicyIndexTemp(:,:,1) .* (mLaborPolicy == 0)  ...
            + mAssetPolicyIndexTemp(:,:,2) .* (mLaborPolicy == 1);

        mValuePath(:,:,iPeriods) = mValueNew;
        mAssetPolicyIndexPath(:,:,iPeriods) =mAssetPolicyIndex;
        mAssetPolicyPath(:,:,iPeriods) =mAssetPolicy;
        mConsumptionPolicyPath(:,:,iPeriods) =mConsumptionPolicy;
        mLaborPolicyPath(:,:,iPeriods) =mLaborPolicy;

    end
    % end
    toc
    fprintf('Iteration finished.\n'); 

    %% FIND TRANSITION CAPITAL & EFFECTIVE LABOR SUPPLIES FORWARD

    mDistributionPath=zeros(nAssets*nGridShocks,nPeriods);
    mDistributionPath(:,1)=vDistributionInitial;
    vAssetPath=zeros(nPeriods,1);
    vAssetPath(1)=capitalSupplyInitial;
    vLaborPath=zeros(nPeriods,1);


    % Find the transition matrices 
    for iPeriods=2:nPeriods

        rows = zeros(1,nGridShocks*nAssets*nGridShocks);
        columns = zeros(1,nGridShocks*nAssets*nGridShocks);
        probs = zeros(1,nGridShocks*nAssets*nGridShocks);
        index = 1;

        for i=1:nAssets % assets
            for j=1:nGridShocks % shock today
                kpIndex = mAssetPolicyIndexPath(i,j,iPeriods-1);
                start = kpIndex;
                for l=1:nGridShocks % shock tomorrow
                    rows(index) = nAssets*(j-1)+i;
                    columns(index) = start+nAssets*(l-1);
                    probs(index) = mTransition(j,l);
                    index=index+1;
                end
            end
        end

        ultimateTransition = sparse(rows,columns,probs,nGridShocks*nAssets,nGridShocks*nAssets);    
        mDistributionPath(:,iPeriods)=(ultimateTransition')*mDistributionPath(:,iPeriods-1);
    end   

    % calculate capital supply path
    for iPeriods=2:nPeriods
      % the capital supply in each period is determined by previous period
      % distribution and policy functions!

        vDistributionAtPrevTime = mDistributionPath(:,iPeriods-1); 
        mAssetPolicyAtPrevTime=mAssetPolicyPath(:,:,iPeriods-1);
        mAssetPolicyAtPrevTimeVectorized=mAssetPolicyAtPrevTime(:);
        vAssetPath(iPeriods)=vDistributionAtPrevTime'*mAssetPolicyAtPrevTimeVectorized;

    end

    % calculate effective labor supply path
    for iPeriods=1:nPeriods
        vDistributionToday = mDistributionPath(:,iPeriods);
        mLaborPolicyToday = mLaborPolicyPath(:,:,iPeriods);
        mLaborPolicyTodayVectorized = mLaborPolicyToday(:);
        vLaborPath(iPeriods)=vDistributionToday' * mLaborPolicyTodayVectorized;
    end

    vKOverLSupplyPath =  vAssetPath./vLaborPath;

    probDistance = norm(mDistributionPath(:,nPeriods)-vDistributionFinal);


    distanceOffStability = max(abs(vKOverLSupplyPath-vKOverLPathNew));
    vDistanceOffStability(iteration+1) = distanceOffStability;
    
    if iteration~=1 % error analysis - as it gets closer be more greedy with g
            % update kOverL path only in the central part

        if distanceOffStability<vDistanceOffStability(iteration)
            vKOverLPathNew(2:nPeriods-1) = 0.5*vKOverLSupplyPath(2:nPeriods-1) + 0.5*vKOverLPathNew(2:nPeriods-1); % conservative update  
        else
            vKOverLPathNew(2:nPeriods-1) = 0.1*vKOverLSupplyPath(2:nPeriods-1) + 0.9*vKOverLPathNew(2:nPeriods-1); % conservative update  
        end
            vKOverLPathNew(vKOverLPathNew>kOverLInitial)=kOverLInitial;     % capital won't be higher than steady state in this particular experiment
    end


    fprintf(' iteration = %d, Stabilization Check = %2.8f, Smoothness Check = %2.8f\n', iteration, distanceOffStability,probDistance)
    iteration=iteration+1;

end



% save withoutUBI_transition.mat


%% 2.  Transition with UBI
% Repeat the same question, but now in the economy with a UBI from above.
% Comment on the differences between the results.

% load withUBI.mat
% This is computationally more complicated than without UBI since now the
% goverment has to make budget balance.