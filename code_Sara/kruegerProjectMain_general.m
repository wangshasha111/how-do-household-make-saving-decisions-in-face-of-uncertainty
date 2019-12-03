% Sara Casella
% Economics 712
% Spring 2019
% Project Dirk Krueger
% General equilibrium: wealth market clearing

clc
clear
close all

%% PART II

%% Question 2.2 
% Reproduce table II from Aiyagari (1994)

% Parameters

% Calibration
N = 11;                   % number of point discretized income process
vDelta = [0,0.3,0.6,0.9]; % income process persistence
sigma = 0.2;             % income process variance
vGamma = [1,3,5];        % 1/EIS
chi = 0;                 % borrowing constraint
alpha = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
tau = 0;                 % tax rate
lambda = 0;              % tax system progressivity
TFP = 1;                 % total factor productivity

parameters.sigma = sigma;
parameters.chi   = chi;
parameters.depreciation  = depreciation; 
parameters.alpha = alpha;
parameters.tau = tau;
parameters.lambda = lambda;

% Computational
nPtsRound1 = 200;
nPtsRound2 = 500;
nPtsRound3 = 3000;
upperBound = 250; % a part from the first case this works well.

% Asset limits
% AssetLimits are deflated by r and turned to the negative
% in partialEqBellman
% Notice pretty small difference for different deltas because of adjustment
% in variance.
AssetLimit=zeros(length(vDelta),1);
for i=1:length(vDelta)
     delta = vDelta(i);
    sigmaEps=sigma*(sqrt(1-delta^2));
    [vIncomeGen, ~]=rouwenhorst(delta,sigmaEps,N);
    vIncomeGen = exp(vIncomeGen)';
    AssetLimit(i)=min(vIncomeGen);
end


% Call function that calculates general equilibrium
% For all the different parameters combinations

general = 1;  % general equilibrium indicator
labourEq = 1; % total labor is 1

options = optimset('Display','Iter','TolX',1e-07);  
equilibriumRate = zeros(length(vDelta),length(vGamma));
difference = zeros(length(vDelta),length(vGamma));
% take initial guesses from Aiyagari (1994) table II
mInitialGuess = [4.15, 4.1365, 4.0912, 3.9305;
                 4.1456, 4.0432, 3.8767, 3.2903;
                 4.0858, 3.9054, 3.5857, 2.4260];
% mInitialGuess = [4.0649, 3.9554, 3.7567, 3.3054;
%                  3.7816, 3.4188, 2.7835, 1.2894;
%                  3.4177, 2.8032, 1.8070, -0.3456];
             
mInitialGuess = mInitialGuess/100;
                
for i=1:length(vDelta)
     for j=1:length(vGamma)
         
        parameters.delta = vDelta(i);
        parameters.gamma = vGamma(j);
        guess = mInitialGuess(j,i);
   
        myfn = @(rate) generalEquilibrium(rate,parameters,TFP,N,nPtsRound1,nPtsRound2,nPtsRound3,...
                                        upperBound,AssetLimit(i),general,labourEq);
        tic                             
        [equilibriumRate(i,j),difference(i,j)] = fzero(myfn,guess,options);
        toc

    end
end

 save('AiyagariTable')
 
 %% Question 2.3
% Tax systems comparisons

clear clc
options = optimset('Display','Iter','TolX',1e-07);  
general = 1;  % general equilibrium indicator
labourEq = 1; % total labor is 1

% Parameters

% Calibration
N = 11;                  % number of point discretized income process
delta = 0.6;             % income process persistence
sigma = 0.2;             % income process variance
gamma = 1;               % 1/EIS
chi = 0;                 % borrowing constraint
alpha = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
TFP = 1;                 % total factor productivity

parameters.delta = delta;
parameters.gamma = gamma;
parameters.sigma = sigma;
parameters.chi   = chi;
parameters.depreciation  = depreciation; 
parameters.alpha  = alpha;

% Computational
nPtsRound1=200;
nPtsRound2=500;
nPtsRound3=3000;
upperBound=160;

% Asset limits
% AssetLimits are deflated by r and turned to the negative
% in partialEqBellman
sigmaEps=sigma*(sqrt(1-delta^2));
[vIncomeGen, mTransitionGen]=rouwenhorst(delta,sigmaEps,N);
vIncomeGen = exp(vIncomeGen)';
AssetLimit=min(vIncomeGen);

% 1) First tax system
tau = 0.18;              % tax rate
lambda = 0.11;           % tax system progressivity

parameters.tau    = tau;
parameters.lambda = lambda;

% Find equilibrium interest rate
myfn = @(rate) generalEquilibrium(rate,parameters,TFP,N,nPtsRound1,nPtsRound2,nPtsRound3,...
                                        upperBound,AssetLimit,general,labourEq);
                                    
guess = 0.0414;

tic                             
eqRateTax = fzero(myfn,guess,options);
toc

% Find total tax renevue
capitalDemand=(alpha/(eqRateTax+depreciation))^(1/(1-alpha))*labourEq; % eq capital
w = (1-alpha)*(capitalDemand^alpha); % eq wage
[vec,~] = eigs(mTransitionGen,1);
vec = abs(vec);
vInvTrans = vec(:)/sum(vec(:)); % invariant distribution
% tax rebated: average tax levied
taxRebated = (w*vIncomeGen-(1-tau)*(w*vIncomeGen).^(1-lambda))'*vInvTrans;

% 2) Second tax system
% Choose τ in such a way that tax revenues (and thus lump-sum rebates to 
% households) are the same under the progressive and the flat tax system.

lambda = 0;           % tax system is not progressive
parameters.lambda = lambda;

myfntax = @(tau) nonProgTaxEquilibrium(tau,parameters,TFP,N,nPtsRound1,nPtsRound2,nPtsRound3,...
    upperBound,AssetLimit,labourEq,eqRateTax,mTransitionGen,vIncomeGen)-taxRebated;

tic                             
optTaxRate = fzero(myfntax,0.21,options);
toc


%% PART III
clc

%optTaxRate = 0.18;
options = optimset('Display','Iter','TolX',1e-07);  
general = 1;

%% Question 3.1
% Transition to a great recession with flat taxes
% Experiment: TFP falls 10% for 10 periods and then goes back to 1.

% Calibration
N = 11;                  % number of point discretized income process
delta = 0.6;             % income process persistence
sigma = 0.2;             % income process variance
gamma = 1;               % 1/EIS
chi = 0;                 % borrowing constraint
alpha = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
tau = optTaxRate;        % tax rate
lambda = 0;              % no progressivity
labourEq = 1;

parameters.delta = delta;
parameters.gamma = gamma;
parameters.sigma = sigma;
parameters.chi   = chi;
parameters.depreciation  = depreciation; 
parameters.alpha  = alpha;
parameters.tau    = tau;
parameters.lambda = lambda;

% Computational
nPtsRound1=200;
nPtsRound2=500;
nPtsRound3=5000;
upperBound=140;
AssetLimit = 0;

% 1) Find steady state equilibrium
TFP = 1; % total factor productivity

% Find steady state interest rate and assets
myfn = @(rate) generalEquilibrium(rate,parameters,TFP,N,nPtsRound1,nPtsRound2,nPtsRound3,...
                                        upperBound,AssetLimit,general,labourEq);
                                 
guess = 0.041266;
tic                         
eqRateSteady = fzero(myfn,guess,options);
toc

[Phi0,Value0,Cons0,Asset0,index0,K0]=equilibriumResults(eqRateSteady,parameters,TFP,...
                                        N,nPtsRound1,nPtsRound2,nPtsRound3,...
                                        upperBound,AssetLimit,general,labourEq);

save('initialDistrData')

 %% Transition
 
load('initialDistrData')

% 2) Guess a sequence of aggregate capital
timeTransition = 200;
%initialK = [linspace(K0,4.7,11),linspace(4.7, K0, timeTransition-11)];
initialK = linspace(K0,K0,timeTransition);

initialPhi         = Phi0;
finalPhi           = Phi0;
initialValue       = Value0;
finalValueFunction = Value0;
initialConsumption = Cons0;
finalConsumption   = Cons0;
initialAsset       = Asset0;
finalAsset         = Asset0;
initialIndex       = index0;
finalIndex         = index0;

% 3) Given initial guess, calculate interest rates and
% find new value functions backwards. Given initial distribution,
% Iterate forward to find new guess. Use relaxation.

nPts = nPtsRound3;
distanceOffStability = 100;
tol = 1e-03;
iteration = 0;
itMax = 100;

K = initialK;
T = length(K);

while distanceOffStability >= tol && iteration <= itMax

   [distanceOffStability, probDistance, mPhi,Asset,mValueFunctionT,...
                            mPolFuncAssetT,mPolFuncConsumptionT,mPolFuncIndexT]=...
                            transitionSolution(N,parameters,...
                          0.0413,nPts, T,...
                          finalValueFunction,initialPhi,finalPhi, K,...
                          upperBound, AssetLimit,finalIndex,...
                          finalAsset,finalConsumption,initialIndex,...
                          initialAsset,initialConsumption,initialValue,K0);
    
    % update capital path only in the central part
    K(2:T-1)=0.1*Asset(2:T-1)+0.9*K(2:T-1); % conservative update  
    K(K>K0)=K0;     % capital won't be higher than steady state in this particular experiment
    
    iteration=iteration+1;

    fprintf(' iteration = %d, Stabilization Check = %2.8f\n', iteration, distanceOffStability)

end
    

distributionConvergence=probDistance;

if distributionConvergence >= tol
    disp 'yes, you must change timeTransition'
else
    disp 'great'
end

% Construct productivity, output and consumption series
TFPpath = ones(1,T);
TFPpath(2:11) = 0.9;
capitalImpulse = K; 
outputImpulse = TFPpath.*capitalImpulse.^alpha;
consumptionImpulse = outputImpulse - [capitalImpulse(2:end),capitalImpulse(end)]+...
    (1-depreciation)*capitalImpulse;

% Change to percentage deviations from ss
capitalImpulsePD = zeros(1,T);
outputImpulsePD   = zeros(1,T);
consumptionImpulsePD = zeros(1,T);
    
for t = 1:T
    capitalImpulsePD(t) = (log(capitalImpulse(t))...
        -log(capitalImpulse(1)))*100;
    outputImpulsePD(t) = (log(outputImpulse(t))...
        -log(outputImpulse(1)))*100;
    consumptionImpulsePD(t) = (log(consumptionImpulse(t))...
        -log(consumptionImpulse(1)))*100;
end

figure
subplot(2,2,1)
plot(TFPpath,'b','LineWidth',2)
xlabel('time')
ylabel('level')
title('TFP (level)')
subplot(2,2,2)
plot(capitalImpulsePD,'b','LineWidth',2)
hold on
plot(zeros(size(capitalImpulsePD)),'--k')
title('Capital')
xlabel('time')
ylabel('% dev ss')
subplot(2,2,3)
plot(outputImpulsePD,'b','LineWidth',2)
hold on
plot(zeros(size(outputImpulsePD)),'--k')
title('Output')
xlabel('time')
ylabel('% dev ss')
subplot(2,2,4)
plot(consumptionImpulsePD,'b','LineWidth',2)
hold on
plot(zeros(size(consumptionImpulsePD)),'--k')
title('Consumption')
xlabel('time')
ylabel('% dev ss')
