%%%%%%%%%%%%%%%%
%%% Function for Value Function Iteration 
% CRRA Utility Function with AR(1) Income Shock, Non Negative Asset 
% With UBI
% Yoshiki Ando 


function [YY,KK,CC,Wage,Interest,EarningDist, IncomeDist, ...
    AssetDist, ConsDist, GiniAsset, LorenzAsset , ...
    GiniIncome, LorenzIncome ,GiniCons, LorenzCons, SocialWelfare ] =...
    MacroAggUBI(Ngrid, Nincome, Eqrr, EqWage, Eqttau, llambda, kkappa,...
  rrho, ssigma,Persist, SigmaY, MaxA  )

%{
% Parameters
Ngrid =300; Nincome=11; Eqrr=0.035; Eqttau=0.15; llambda=0.2; kkappa=1;...
  rrho=0.04; ssigma=1;Persist=0.8; SigmaY=0.2; MaxA=30;
%}
  
  
% Compute Wage and Capital given rr
aalpha = 0.36; depreciate = 0.08;
Wwage  = @(rr) (1-aalpha) * (aalpha / (rr+depreciate) )^(aalpha/(1-aalpha)) ;
% Wage   = Wwage(Eqrr) ;
Wage = EqWage ;
Kcapital = @(rr) (aalpha / (rr+depreciate) )^(1/(1-aalpha)) ;
KK  = Kcapital (Eqrr) ;

% Solve Value Function Iteration Given Parameters
[mValue, mAssetLocation,mAssetChoice, mLaborChoice, mConsChoice] =...
    MultigridUBI(Ngrid, Nincome, Eqrr, Wage, ...
    Eqttau, llambda, kkappa,rrho,  ssigma,Persist, SigmaY, MaxA  );


% Compute Stationary Distribution given Parameters
[vJointStationary, StationaryMat ,vStationaryAsset,mLaborChoice] =...
    DistUBI(Ngrid, Nincome, Eqrr, Eqttau, llambda, kkappa, rrho,...
    ssigma,Persist, SigmaY, MaxA  ) ;

% Compute Aggregate Labor Supply
LL = sum( StationaryMat .* mLaborChoice, 'all' ) ;
YY = KK^aalpha * LL^(1-aalpha)  ;
CC = sum( StationaryMat .* mConsChoice, 'all' ) ;
Interest = Eqrr;

% Earnings  Distribution
% Tauchen
[vState, ~, vStatIncome] = myTauchenAR1(Persist, 0, SigmaY, Nincome);

% Probability for agents who work
vEarningsTemp    = (1-Eqttau) .* Wage .* exp(vState) ./ (exp(vState) * vStatIncome);
vEarningDistTemp =  sum(StationaryMat.* mLaborChoice  ,1) ;

% Include the probability of agents who don't work (earnings = 0)
vEarnings     = [0,vEarningsTemp] ;
vEarningDist  = [ sum(StationaryMat.* (1-mLaborChoice) ,'all'), vEarningDistTemp] ;
EarningDist   = [vEarnings' , vEarningDist' ] ;

% Preparation for titles in Figures
 if llambda ==0
    TitleUBI = 'without UBI' ;
 else
    TitleUBI = 'under UBI' ; 
 end

title1 = strcat('Distribution of Earnings' , {' '}, TitleUBI) ;
 
figure;
pl1 = bar(EarningDist(:,1), EarningDist(:,2) ) ;
tit = title(title1, 'FontSize' , 18) ;


% Income Distribution
vCurrentAsset    = linspace(0,MaxA,Ngrid) ;
mIncome          = repmat( vEarningsTemp, [Ngrid,1] ) .* mLaborChoice +...
    Eqrr.* repmat( vCurrentAsset', [1,Nincome] ) ;
vectorizeIncome     = reshape(mIncome,numel(mIncome),1) ;
vectorizeIncomeProb = reshape(StationaryMat,numel(StationaryMat),1);

% Round and combine the element with the same value
vectorizeIncome2 = round(vectorizeIncome,2) ;
[a1,~,c1] = unique(vectorizeIncome2);
A1 = accumarray(c1,vectorizeIncomeProb);
IncomeDist = [a1,A1]   ;
title2 = strcat('Distribution of Income' , {' '}, TitleUBI) ;
figure;
pl2 = bar(IncomeDist(:,1), IncomeDist(:,2), 'BarWidth', 0.85) ;
ti2=title(title2, 'FontSize',18);

% Asset Distribution
vStationaryAsset  = sum(StationaryMat,2) ;
AssetDist   = [vCurrentAsset', vStationaryAsset] ;
title13 = strcat('Distribution of Asset' , {' '}, TitleUBI) ;
figure;
pl3 = bar(AssetDist(:,1), AssetDist(:,2)) ;
ti3=title(title13 , 'FontSize', 18);

% Consumption Distribution
vectorizeCons     = reshape(mConsChoice,numel(mConsChoice),1) ;
vectorizeConsProb = reshape(StationaryMat,numel(StationaryMat),1);

% Round and combine the element with the same value
vectorizeIncome2 = round(vectorizeCons,2) ;
[a2,~,c2] = unique(vectorizeIncome2);
A2 = accumarray(c2,vectorizeConsProb);
ConsDist = [a2,A2]   ;

title4 = strcat('Distribution of Consumption' , {' '}, TitleUBI) ;
figure;
pl4 = bar(ConsDist(:,1), ConsDist(:,2), 'BarWidth', 0.85) ;
ti4=title(title4, 'FontSize', 18);

% Lorenz Curve and Gini Coefficient for Asset
AssetDist          = [vCurrentAsset', vStationaryAsset] ;
vAssetProbCum      = round( cumsum(vStationaryAsset),3) ;
vAssetTimesProb    = vCurrentAsset' .* vStationaryAsset ;
vAssetTimesProbCum = cumsum(vAssetTimesProb) / sum(vAssetTimesProb) ;
LorenzAsset        = [vAssetProbCum, vAssetTimesProbCum] ;

GiniAsset = (0.5 - trapz(LorenzAsset(:,1), LorenzAsset(:,2)) )*2 ;

% Lorenz Curve and Gini Coefficient for Income
vIncomeProbCum      = round( cumsum(IncomeDist(:,2)),3) ;
vIncomeTimesProb    = IncomeDist(:,1) .* IncomeDist(:,2) ;
vIncomeTimesProbCum = cumsum(vIncomeTimesProb) / sum(vIncomeTimesProb) ;
LorenzIncome        = [ [0;vIncomeProbCum], [0;vIncomeTimesProbCum] ] ;

GiniIncome = (0.5 - trapz(LorenzIncome(:,1), LorenzIncome(:,2)) )*2 ;

% Lorenz Curve and Gini Coefficient for Consumption
vConsProbCum      = round( cumsum(ConsDist(:,2)),3) ;
vConsTimesProb    = ConsDist(:,1) .* ConsDist(:,2) ;
vConsTimesProbCum = cumsum(vConsTimesProb) / sum(vConsTimesProb) ;
LorenzCons        = [vConsProbCum, vConsTimesProbCum] ;

GiniCons = (0.5 - trapz(LorenzCons(:,1), LorenzCons(:,2)) )*2 ;



% Figure for three Lorenz Curves
title15 = strcat('Lorenz Curve' , {' '}, TitleUBI) ;
figure;
Uniline = linspace(0,1,100) ;
pl = plot( LorenzAsset(:,1), LorenzAsset(:,2),LorenzIncome(:,1),...
    LorenzIncome(:,2), LorenzCons(:,1), LorenzCons(:,2),Uniline,Uniline);
tit=title(title15, 'FOntSize',18);
set(pl,{'LineStyle'},{'-';'-.'; '--';':'});
set(pl,'Linewidth',2);
le=legend({'Asset','Income','Consumption','45 degree line'},...
    'Location','northwest','FontSize',12);

% SocialWelfare
SocialWelfare = sum( StationaryMat .* mValue, 'all' ) ;


end