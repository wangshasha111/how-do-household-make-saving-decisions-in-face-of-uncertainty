function[vGridAsset,vIncome,mTransition,mValueFunction,mPolFuncAsset,...
        mPolFuncConsumption,mPolFuncIndex,nGridAsset]=partialEqBellman(N,capitalDemand,...
                              parameters,interestRate,nPtsRound1,nPtsRound2,nPtsRound3,...
                              upperBound, AssetLimit, general)

% OUTPUT
%--------------------------------------------------------------------------
% vGridAsset: bond grid
%--------------------------------------------------------------------------
% vIncome: discretized income shock grid
%--------------------------------------------------------------------------
% mTransition: transition matrix
%--------------------------------------------------------------------------
% mPolFuncIndex: a policy function storing the values of the indexes
% corresponding to optimal bond choice, given old bond level and the
% realization of the shock
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% N: number of element used in the Rouwen-Horst discretization of the
% income process
%--------------------------------------------------------------------------
% Capital Demand (USED ONLY FOR GENERAL EQM EXERCISES): uses the given
% parameters and interest rate to derive the corresponding demand for
% capital on the firm side
%--------------------------------------------------------------------------
% parameters is an array that must include the following:

% delta: persistence of the income process
% sigma: variance of the income process
% gamma: elasticity of substitution parameter
% rho: rate of time preference
% chi: tightness of the borrowing constraint
% alpha (USED ONLY FOR GENERAL EQM EXERCISES): capital share
% tau (USED ONLY FOR GENERAL EQM EXERCISES): tax rate
% lambda (USED ONLY FOR GENERAL EQM EXERCISES): tax progressivity
%--------------------------------------------------------------------------
% interestRate: interest rate
%--------------------------------------------------------------------------
% nPtsRound1, nPtsRound2, nPtsRound3: grid points for the 3 multigrid
% step taken to approximate Bellman's solution
%--------------------------------------------------------------------------
% upperBound, AssetLimit: for the bond grid
%--------------------------------------------------------------------------
% general: dichotomous variable which is zero if working in partial
% equilibrium and one otherwise. In partial equilibrium r is an exogenous
% parameter that takes part in saving determination and income is totally
% exogenous; in general equilibrium there is a market wage component and
% supply and demand of capital interact through r.
%--------------------------------------------------------------------------

%% Recover parameters

a = 0.2;    % spacing of the grid (exponential)

delta = parameters.delta; % persistence of the income process
sigma = parameters.sigma; % variance of the income process
gamma = parameters.gamma; % inverse elasticity of substitution
chi   = parameters.chi;   % tightness of the borrowing constraint

r     = interestRate;     % interest rate

if general == 0
    rho  = parameters.rho;   % rate of time preference
    beta = 1/(1+rho);
else
    rho = 0.0417;
    beta = 1/(1+rho); 
    %alpha = parameters.alpha;   % (USED ONLY FOR GENERAL EQM EXERCISES) capital share
    w = parameters.wage;
    tau = parameters.tau;       % tax rate 
    lambda = parameters.lambda; % tax progressivity
end
%% Income Process

  if general == 0
      sigma=sigma*(sqrt(1-delta^2));
      [vIncome, mTransition]=rouwenhorst(delta,sigma,N);
      vIncome = exp(vIncome);
      vIncome=vIncome';
  else
      
      sigma=sigma*(sqrt(1-delta^2));
      [vIncome,mTransition]=rouwenhorst(delta,sigma,N);
      vIncome=vIncome';
      vIncome = exp(vIncome);
      [vec,~] = eigs(mTransition,1);
      vec = abs(vec);
      vInvTrans = vec(:)/sum(vec(:));     
      % tax rebated: average tax levied
      taxRebated = (w*vIncome-(1-tau)*(w*vIncome).^(1-lambda))'*vInvTrans;

  end
  
  % If we use rouwen the probability matrix that is outputted must not be
  % transposed in the value function iteration part, and the sigma is the
  % sigma of the income process
  
  % since the markov chain is on the log of income, log(y')=delta*log(y)+
  %                                                        (1-delta^2)eps'
  



%% Grid Generation

minK=-chi*AssetLimit/r;
maxK=upperBound;

%vGridAsset = linspace(minK, maxK, nPtsRound1);
temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nPtsRound1));
temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits

nGridAsset = length(vGridAsset);
[~, index]=min(abs(vGridAsset));
vGridAsset(index)=0;

nGridStoch=length(vIncome);

%% Required Matrices and Vectors

% Initial Value Functions

mValueFunction    = zeros(nGridAsset,nGridStoch); % using 0 as first Guess (average utility)
mValueFunctionNew = zeros(nGridAsset,nGridStoch);
mPolFuncAsset     = zeros(nGridAsset,nGridStoch);
mPolFuncConsumption = zeros(nGridAsset,nGridStoch);
mPolFuncIndex = zeros(nGridAsset,nGridStoch);

%% Value function iteration

if general==0

    maxDifference = 10.0;
    tolerance = 10^-8;
    iteration = 0;

    while (maxDifference>tolerance)  

         expectedValueFunction = mValueFunction*(mTransition');


      for iIncome=1:nGridStoch


        z=vIncome(iIncome,1); 
        y=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset

                Asset=vGridAsset(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset
                    Asset1=vGridAsset(iAssetTplus1);
                    Consumption=(y+(1+r)*Asset-Asset1);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption(iAsset,iIncome) = (z+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
        mValueFunction = mValueFunctionNew;

        iteration = iteration+1;


    end

else
    
    maxDifference = 10.0;
    tolerance = 10^-8;

    while (maxDifference>tolerance)  

         expectedValueFunction = mValueFunction*(mTransition');

      for iIncome=1:nGridStoch

        z=vIncome(iIncome,1);
        l=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset

                Asset=vGridAsset(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset
                    Asset1=vGridAsset(iAssetTplus1);
                    Consumption=((1-tau)*(w*l)^(1-lambda)+(1+r)*Asset-Asset1+taxRebated);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption(iAsset,iIncome) = (w*l+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
        mValueFunction = mValueFunctionNew;


    end
end

%% Multigrid Round 2

%vGridAsset2 = linspace(minK, maxK, nPtsRound2);
temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nPtsRound2));
temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
vGridAsset2 = temp_01*(maxK-minK) + minK; % re-scale to be between limits

nGridAsset2 = length(vGridAsset2);
mValueFunction2      = zeros(nGridAsset2,nGridStoch);
mValueFunctionNew2   = zeros(nGridAsset2,nGridStoch);
mPolFuncAsset2       = zeros(nGridAsset2,nGridStoch);
mPolFuncConsumption2 = zeros(nGridAsset2,nGridStoch);
mPolFuncIndex2 = zeros(nGridAsset2,nGridStoch);

for i=1:nGridStoch
    mValueFunction2(:,i)=interpn(vGridAsset, mValueFunction(:,i), vGridAsset2);
end

if general==0

    maxDifference = 10.0;
    tolerance = 10^-8;
    iteration = 0;

    while (maxDifference>tolerance)  

         expectedValueFunction2 = mValueFunction2*(mTransition');


      for iIncome=1:nGridStoch


        z=vIncome(iIncome,1); 
        y=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset2
                Asset=vGridAsset2(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset2  

                    Asset1=vGridAsset2(iAssetTplus1);
                    Consumption=(y+(1+r)*Asset-Asset1);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction2(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction2(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset2(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex2(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew2(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption2(iAsset,iIncome) = (z+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew2-mValueFunction2)));
        mValueFunction2 = mValueFunctionNew2;

        iteration = iteration+1;


    end

else
    
    maxDifference = 10.0;
    tolerance = 10^-8;

    while (maxDifference>tolerance)  

         expectedValueFunction2 = mValueFunction2*(mTransition');


      for iIncome=1:nGridStoch

        z=vIncome(iIncome,1);
        l=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset2

                Asset=vGridAsset2(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset2
                    Asset1=vGridAsset2(iAssetTplus1);
                    Consumption=((1-tau)*(w*l)^(1-lambda)+(1+r)*Asset-Asset1+taxRebated);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction2(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction2(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset2(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex2(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew2(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption2(iAsset,iIncome) = (w*l+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew2-mValueFunction2)));
        mValueFunction2 = mValueFunctionNew2;


    end
end
%% Multigrid Round 3

% vGridAsset3 = linspace(minK, maxK, nPtsRound3);
temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nPtsRound3));
temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
vGridAsset3 = temp_01*(maxK-minK) + minK; % re-scale to be between limits

nGridAsset3 = length(vGridAsset3);
mValueFunction3 = zeros(nGridAsset3,nGridStoch);
mValueFunctionNew3 = zeros(nGridAsset3,nGridStoch);
mPolFuncAsset3   = zeros(nGridAsset3,nGridStoch);
mPolFuncConsumption3   = zeros(nGridAsset3,nGridStoch);
mPolFuncIndex3 = zeros(nGridAsset3,nGridStoch);


for i=1:nGridStoch
    mValueFunction3(:,i)=interpn(vGridAsset2, mValueFunction2(:,i), vGridAsset3);
end

if general==0

    maxDifference = 10.0;
    tolerance = 10^-8;
    iteration = 0;

    while (maxDifference>tolerance)  

         expectedValueFunction3 = mValueFunction3*(mTransition');

      for iIncome=1:nGridStoch


        z=vIncome(iIncome,1); 
        y=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset3

                Asset=vGridAsset3(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset3
                    Asset1=vGridAsset3(iAssetTplus1);
                    Consumption=(y+(1+r)*Asset-Asset1);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction3(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction3(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset3(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex3(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew3(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption3(iAsset,iIncome) = (z+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew3-mValueFunction3)));
        mValueFunction3 = mValueFunctionNew3;

        iteration = iteration+1;

    end

else

    maxDifference = 10.0;
    tolerance = 10^-8;

    while (maxDifference>tolerance)  

         expectedValueFunction3 = mValueFunction3*(mTransition');

      for iIncome=1:nGridStoch

        z=vIncome(iIncome,1);
        l=z;

        gridAssetTplus1 = 1; %including non-zero asset holding constraint

        for iAsset=1:nGridAsset3

                Asset=vGridAsset3(iAsset);
                valueHighSoFar = -Inf;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset3
                    Asset1=vGridAsset3(iAssetTplus1);
                    Consumption=((1-tau)*(w*l)^(1-lambda)+(1+r)*Asset-Asset1+taxRebated);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(max(Consumption,eps)))+...
                                       beta*expectedValueFunction3(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((max(Consumption,eps)^(1-gamma)-1)/(1-gamma)))+...
                                       beta*expectedValueFunction3(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAsset3(iAsset,iIncome) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        mPolFuncIndex3(iAsset,iIncome) = gridAssetTplus1;
                    else
                        break; % We break when we have achieved the max
                    end    

                end


                mValueFunctionNew3(iAsset,iIncome) = valueHighSoFar;
                mPolFuncConsumption3(iAsset,iIncome) = (w*l+(1+r)*Asset-Asset1);


        end

       end

        maxDifference = max(max(abs(mValueFunctionNew3-mValueFunction3)));
        mValueFunction3 = mValueFunctionNew3;


    end

end

mValueFunction = mValueFunction3;
mPolFuncAsset  = mPolFuncAsset3;
mPolFuncConsumption = mPolFuncConsumption3;  
mPolFuncIndex = mPolFuncIndex3;
vGridAsset = vGridAsset3;
nGridAsset    = nGridAsset3;

end
