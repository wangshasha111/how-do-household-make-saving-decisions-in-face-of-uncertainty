function[distanceOffStability, probDistance, mPhi,Asset,mValueFunctionT,...
           mPolFuncAssetT,mPolFuncConsumptionT,mPolFuncIndexT]=transitionSolution(N,parameters,...
                          eqRateSteady,nPts, T,...
                          finalValueFunction,initialPhi,finalPhi, K,...
                          upperBound, AssetLimit,finalIndex,...
                          finalAsset,finalConsumption,initialIndex,...
                          initialAsset,initialConsumption,initialValue,initialK)

%% Parameters   

rho = 0.0417;
beta = 1/(1+rho);

a = 0.2;    % spacing of the grid (exponential)

delta = parameters.delta; % persistence of the income process
sigma = parameters.sigma; % variance of the income process
gamma = parameters.gamma; % inverse elasticity of substitution
chi   = parameters.chi;   % tightness of the borrowing constraint
tau    = parameters.tau;    % tax rate 
lambda = parameters.lambda; % tax progressivity
alpha = parameters.alpha;
depreciation = parameters.depreciation;

%% Income Process

sigma=sigma*(sqrt(1-delta^2));
[vInT,mTransitionT]=rouwenhorst(delta,sigma,N);
vInT=vInT';
vInT = exp(vInT);
[vec,~] = eigs(mTransitionT,1);
vec = abs(vec);
vInvTrans = vec(:)/sum(vec(:));  

%% Grid Generation

r = eqRateSteady;
%making an approximating assumption and taking this constraint as fixed
%during transition, at the final interest level.
minK=-chi*AssetLimit/r;
maxK=upperBound;

%vGridAsset = linspace(minK, maxK, nPtsRound1);
temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nPts));
temp_01 = temp/max(temp) - min(temp)/max(temp); % re-scale to be between 0 and 1
vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
nGridAsset = length(vGridAsset);

%[~, index]=min(abs(vGridAsset));
%vGridAsset(index)=0;

nGridStoch=length(vInT);

vTime=(1:1:T);
nTime=length(vTime);


%% Required Matrices and Vectors

% Initial Value FunctionTs

mValueFunctionT = zeros(nGridAsset,nGridStoch,nTime);
mValueFunctionT(:,:,T) = finalValueFunction;
mPolFuncAssetT   = zeros(nGridAsset,nGridStoch,nTime);
mPolFuncAssetT(:,:,T) = finalAsset;
mPolFuncConsumptionT   = zeros(nGridAsset,nGridStoch,nTime);
mPolFuncConsumptionT(:,:,T) = finalConsumption;
mPolFuncIndexT = zeros(nGridAsset,nGridStoch,nTime);
mPolFuncIndexT(:,:,T) = finalIndex;

%% FIND VALUE FUNCTION ITERATING BACKWARDS

for time=nTime-1:-1:2
     
     capitalDemand=K(time);
      if time > 1 && time <= 11
          TFP = 0.9;
      else
          TFP = 1;
      end
%    TFP = 0.98;
     wTime=(1-alpha)*TFP*(capitalDemand^alpha);
     rTime=alpha*TFP*(capitalDemand^(alpha-1))-depreciation;
     
     taxRebated = (wTime*vInT-(1-tau)*(wTime*vInT).^(1-lambda))'*vInvTrans;
     expectedValueFunctionT = mValueFunctionT(:,:,time+1)*mTransitionT';
 
      for iIncome=1:nGridStoch

        l=vInT(iIncome,1);

        gridAssetTplus1 = 1;

        for iAsset=1:nGridAsset

                Asset=vGridAsset(iAsset);
                 valueHighSoFar = -1000.0;

                for iAssetTplus1 = gridAssetTplus1:nGridAsset
                    Asset1=vGridAsset(iAssetTplus1);
                    Consumption=((1-tau)*(wTime*l)^(1-lambda)+(1+rTime)*Asset-Asset1+taxRebated);

                    if gamma==1 %log utility case
                      valueProvisional = (1-beta)*(log(Consumption))+...
                                       beta*expectedValueFunctionT(iAssetTplus1,iIncome);
                    else
                      valueProvisional = (1-beta)*(((Consumption)^(1-gamma)-1)/(1-gamma))+...
                                       beta*expectedValueFunctionT(iAssetTplus1,iIncome);
                    end

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        mPolFuncAssetT(iAsset,iIncome,time) = Asset1;
                        gridAssetTplus1 = iAssetTplus1;
                        assetIndex= iAssetTplus1;

                    else
                        break; % We break when we have achieved the max
                    end    

                end

                mPolFuncIndexT(iAsset,iIncome,time) = assetIndex;      
                mValueFunctionT(iAsset,iIncome,time) = valueHighSoFar;
                mPolFuncConsumptionT(iAsset,iIncome,time) = Consumption; 

         end

      end
       
end

% Forcing the policy index to conform with initial and final SS
mPolFuncIndexT(:,:,1)=initialIndex;
mPolFuncAssetT(:,:,1)=initialAsset;
mPolFuncConsumptionT(:,:,1)=initialConsumption;
mValueFunctionT(:,:,1)=initialValue;
mPolFuncIndexT(:,:,T-1) = finalIndex;
mPolFuncAssetT(:,:,T-1)=finalAsset;


%% FIND TRANSITION CAPITAL SUPPLIES FORWARD

mPhi=zeros(nGridAsset*N,T);
mPhi(:,1)=initialPhi;
Asset=zeros(1,T);
Asset(1,1)=initialK;

% Find the transition matrices 
for transPeriod=2:T
    
%      ultimateTransition = zeros(nGridAsset*N,nGridAsset*N);
%      for j = 1:N % income today
%          for k = 1:N % income tomorrow
%             for i = 1:nGridAsset % assets tomorrow
%                     ultimateTransition((j-1)*nGridAsset+i,(k-1)*nGridAsset+...
%                         mPolFuncIndexT(i,j,transPeriod-1)) = mTransitionT(j,k);
%             end
%          end
%      end

%      rows = zeros(1,N*nGridAsset*N);
%      columns = zeros(1,N*nGridAsset*N);
%      probs = zeros(1,N*nGridAsset*N);
%      index = 1;
%      mTransitionT = mTransitionT';
% 
%      for i=1:nGridAsset
%         for j=1:N
%           kpIndex = mPolFuncIndexT(i,j,transPeriod-1);
%           start = N*(kpIndex-1);
%            for l=1:N
%               rows(index) = N*(i-1)+j;
%               columns(index) = start+l;
%               probs(index) = mTransitionT(j,l);
%               index=index+1;
%            end
%         end
%      end
% 
%    ultimateTransition = sparse(rows,columns,probs,N*nGridAsset,N*nGridAsset);

rows = zeros(1,N*nGridAsset*N);
columns = zeros(1,N*nGridAsset*N);
probs = zeros(1,N*nGridAsset*N);
index = 1;
 
for i=1:nGridAsset % assets
 for j=1:N % shock today
   kpIndex = mPolFuncIndexT(i,j,transPeriod-1);
   start = kpIndex;
    for l=1:N % shock tomorrow
       rows(index) = nGridAsset*(j-1)+i;
       columns(index) = start+nGridAsset*(l-1);
       probs(index) = mTransitionT(j,l);
       index=index+1;
    end
 end
end
 
ultimateTransition = sparse(rows,columns,probs,N*nGridAsset,N*nGridAsset);
    
mPhi(:,transPeriod)=(ultimateTransition')*mPhi(:,transPeriod-1);


end   
        
for transPeriod=2:T
  % the capital supply in each period is determined by previous period
  % distribution and policy functions!
    
    mPhiAtPrevTime = mPhi(:,transPeriod-1); 
    mPolFuncAssetAtPrevTime=mPolFuncAssetT(:,:,transPeriod-1);
    mPolFuncAssetVectorized=mPolFuncAssetAtPrevTime(:);
    Asset(1,transPeriod)=mPhiAtPrevTime'*mPolFuncAssetVectorized;

   % Asset(1,transPeriod)=sum(sum(reshape(mPhi(:,transPeriod),[N,nGridAsset]).*mPolFuncAssetT(:,:,transPeriod)'));
    
end
    
%   mPolFuncAssetTCurrent=mPolFuncAssetT(:,:,transPeriod);
%   mPolFuncAssetVectorized=mPolFuncAssetTCurrent(:);
%   Asset(1,transPeriod)=mPhi(:,transPeriod)'*mPolFuncAssetVectorized;

distanceOffStability = max(abs(Asset-K));
  
probDistance = norm(mPhi(:,T)-finalPhi);

end
