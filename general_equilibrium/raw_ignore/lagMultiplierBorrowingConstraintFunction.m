function lagMultiplierBorrowingConstraintInConsumptionUnits = lagMultiplierBorrowingConstraintFunction(mAssetPolicyIndex,mConsumptionPolicy, ggamma, nGridAsset, nGridShocks, vIncomeShocks, mTransition, vGridAsset, bbeta, r)
% stochastic Euler equation
% c^(-ggamma) = bbeta * (1 + r) * E(c' ^(-ggama)) + mmu/(1-bbeta); % Note mmu>=0 is the lagrangian multiplier of the borrowing constraint a'>=o
% stochastic Euler equation normalized to be in the consumption unit
% 1 = c^(ggamma)  * bbeta * (1 + r) * E(c' ^(-ggama)) + c^(ggamma) * mmu/(1-bbeta); 
% 1 - bbeta * (1 + r) * c^(ggamma)  * E(c' ^(-ggama)) == c^(ggamma) * mmu/(1-bbeta); 

    marginalUtilityOfConsumptionToday = mConsumptionPolicy.^(-ggamma);
    
    expectedMarginalUtilityOfConsumptionTomorrow = zeros(nGridAsset,nGridShocks); 
    
    for iShocks = 1 : nGridShocks
%         y=vIncomeShocks(iShocks); 
        for iAsset = 1 : nGridAsset
%             asset=vGridAsset(iAsset);
            iAssetPrime = mAssetPolicyIndex(iAsset,iShocks);
            
            expected = zeros(nGridShocks,1);
            
            for iShocksPrime = 1:nGridShocks
%                 yPrime=vIncomeShocks(iShocksPrime); 

    %             kPrimePrime = mKPolicy(ikPrime,iaPrime);
    %             consumptionPrime_1 = mConsumptionPolicy_1(ikPrime,iaPrime);
    %             consumptionPrime_2 = mConsumptionPolicy_2(ikPrime,iaPrime);
    %             laborPrime_1 = mLaborPolicy_1(ikPrime,iaPrime);
                consumptionPrime = mConsumptionPolicy(iAssetPrime,iShocksPrime);
                marginalUtilityTomorrow = consumptionPrime ^ (-ggamma);
                unexpected = bbeta * (1+r) * marginalUtilityTomorrow;
                expected(iShocksPrime) = unexpected * mTransition(iShocks,iShocksPrime);

            end
            expectedMarginalUtilityOfConsumptionTomorrow(iAsset,iShocks) = sum(expected);
        end % asset
    end % shock
        

% expectedMarginalUtilityOfConsumptionTomorrow = .^(-ggamma) * mTransition';

    lagMultiplierBorrowingConstraintInConsumptionUnits = 1 - expectedMarginalUtilityOfConsumptionTomorrow./marginalUtilityOfConsumptionToday;

end