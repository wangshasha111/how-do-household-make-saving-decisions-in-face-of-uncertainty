function distanceEquilibrium=generalEquilibrium(rate,parameters,TFP,N,nPtsRound1,nPtsRound2,nPtsRound3,...
                                                upperBound,AssetLimit,general,labourEq)

% Recover parameters
alpha = parameters.alpha; % capital share
depreciation = parameters.depreciation; % capital depreciation

% 1) Calculate demand for capital given r
capitalDemand=(alpha*TFP/(rate+depreciation))^(1/(1-alpha))*labourEq;

% 2) Calculate wage given r and demand for capital
wage = (1-alpha)*TFP*(capitalDemand^alpha); % wage
parameters.wage = wage;

% tic
% [~,~,mTransition,~,mPolFuncAsset,~,mPolFuncIndex,nGridAsset]=...
%     partialEqBellmanInterp2(N,parameters,rate,nPtsRound1,nPtsRound2,nPtsRound3,...
%     upperBound, AssetLimit)
% toc
% 3) Solve the household problem for given r and W(r)
[~,~,mTransition,~,mPolFuncAsset,~,mPolFuncIndex,nGridAsset]=partialEqBellman(N,capitalDemand,...
                                            parameters,rate,nPtsRound1,nPtsRound2,nPtsRound3,...
                                            upperBound,AssetLimit,general);


                                        
% % 4) Calculate the total assets in the economy
% ultimateTransition = zeros(nGridAsset*N,nGridAsset*N);
%  for j = 1:N % income today
%      for k = 1:N % income tomorrow
%         for i = 1:nGridAsset % assets tomorrow
%                 ultimateTransition((j-1)*nGridAsset+i,(k-1)*nGridAsset+...
%                     mPolFuncIndex(i,j)) = mTransition(j,k);
%         end
%      end
%  end
% 
%  ultimateTransition = sparse(ultimateTransition);
% 
%      rows = zeros(1,N*nGridAsset*N);
%      columns = zeros(1,N*nGridAsset*N);
%      probs = zeros(1,N*nGridAsset*N);
%      index = 1;
%      mTransition = mTransition';
% % 
%      for i=1:nGridAsset
%         for j=1:N
%           kpIndex = mPolFuncIndex(i,j);
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

rows = zeros(1,N*nGridAsset*N);
columns = zeros(1,N*nGridAsset*N);
probs = zeros(1,N*nGridAsset*N);
index = 1;
 
for i=1:nGridAsset % assets
 for j=1:N % shock today
   kpIndex = mPolFuncIndex(i,j);
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
mProbConverged = vecP/sum(sum(vecP)); 

mPolFuncAssetVectorized=mPolFuncAsset(:);
capitalSupply=mProbConverged'*mPolFuncAssetVectorized;
%capitalSupply = sum(sum(reshape(mProbConverged,[N,nGridAsset]).*mPolFuncAsset'));

 % 5) Check whether total capital = total assets
distanceEquilibrium=capitalDemand-capitalSupply;

end