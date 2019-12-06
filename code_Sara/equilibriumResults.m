function[eqPhi,eqValue,eqCons,eqAsset,eqIndex,eqK]=equilibriumResults(rate,parameters,TFP,...
                                        N,nPtsRound1,...
                                        nPtsRound2,nPtsRound3,...
                                        upperBound,AssetLimit,general,labourEq)
  
% Recover parameters
alpha = parameters.alpha; % capital share
depreciation = parameters.depreciation; % capital depreciation

% Total capital
capitalDemand=(alpha*TFP/(rate+depreciation))^(1/(1-alpha))*labourEq;

% Wage
wage = (1-alpha)*TFP*(capitalDemand^alpha); % wage
parameters.wage = wage;

% Partial equilibrium problem
[vGridAsset,vIncome,mTransition, eqValue, eqAsset,eqCons,eqIndex,nGridAsset]=partialEqBellman(N,capitalDemand,...
                              parameters,rate,nPtsRound1,nPtsRound2,nPtsRound3,...
                              upperBound, AssetLimit, general);


% Calculate the total assets in the economy

%     ultimateTransition = zeros(nGridAsset*N,nGridAsset*N);
%      for j = 1:N % income today
%          for k = 1:N % income tomorrow
%             for i = 1:nGridAsset % assets tomorrow
%                     ultimateTransition((j-1)*nGridAsset+i,(k-1)*nGridAsset+...
%                         eqIndex(i,j)) = mTransition(j,k);
%             end
%          end
%      end

%      rows = zeros(1,N*nGridAsset*N);
%      columns = zeros(1,N*nGridAsset*N);
%      probs = zeros(1,N*nGridAsset*N);
%      index = 1;
%      mTransition = mTransition';
% 
%      for i=1:nGridAsset
%         for j=1:N
%           kpIndex = eqIndex(i,j);
%           start = N*(kpIndex-1);
%            for l=1:N
%               rows(index) = N*(i-1)+j;
%               columns(index) = start+l;
%               probs(index) = mTransition(j,l);
%               index=index+1;
%            end
%         end
%      end
% 
%     ultimateTransition = sparse(rows,columns,probs,N*nGridAsset,N*nGridAsset);
%          
% [vecP,~] = (eigs(ultimateTransition',1));
% vecP = abs(vecP);
% eqPhi = vecP/sum(sum(vecP)); 

 %  mPolFuncAssetVectorized=eqAsset(:);
 %  eqK=eqPhi'*mPolFuncAssetVectorized;
 %eqK = sum(sum(reshape(eqPhi,[N,nGridAsset]).*eqAsset'));
 
rows = zeros(1,N*nGridAsset*N);
columns = zeros(1,N*nGridAsset*N);
probs = zeros(1,N*nGridAsset*N);
index = 1;
 
for i=1:nGridAsset % assets
 for j=1:N % shock today
   kpIndex = eqIndex(i,j);
   start = kpIndex;
    for l=1:N % shock tomorrow
       rows(index) = nGridAsset*(j-1)+i;
       columns(index) = start+nGridAsset*(l-1);
       probs(index) = mTransition(j,l);
       index=index+1;
    end
 end
end
 
ultimateTransition = sparse(rows,columns,probs,N*nGridAsset,N*nGridAsset);

[vecP,~] = (eigs(ultimateTransition',1,'lr','Tolerance',1e-8));
vecP = abs(vecP);
eqPhi = vecP/sum(sum(vecP)); 

mPolFuncAssetVectorized=eqAsset(:);
eqK=eqPhi'*mPolFuncAssetVectorized;
   
                                       
   end