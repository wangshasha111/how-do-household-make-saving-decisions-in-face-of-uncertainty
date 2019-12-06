%% Question 3.3.1
% Evaluate Universal Basic Income

%% Step 1
% Calibrate the economy without the UBI (that is, ¦Ë = ¦Ó = 0) has 80%
% of the population working, that is, find ¦Ê such that in the associated stationary
% equilibrium

%% How I did it:
% Given kkappa, I used the fminbnd to find K/L
% Then in the outer loop, I used the fminbnd to find kkappa to let labor
% participation to be 0.8

clear

%% Grid for labor - discrete labor choice
ifLabor = [0,1];
nGridLabor = length(ifLabor);

% households
rrho = 0.04; % discount rate
bbeta = 1/(1+rrho); % discount factor
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.95; % persistence of the income process
% ssigmaY = [0.2,0.8]; % income process variance
ssigmaY = 0.4; % income process variance

% firms
aalphaK = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
TFP = 1; % total production factor

% income shocks
nGridShocks = 31;
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
upperBound = 10 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?

% Grid of assets
a = 0.2;    % spacing of the grid (exponential)
nAssets = 100;

% Optimization
% options = optimset('Display','Iter','TolX',1e-07);  
options = optimset('Display', 'off','TolX',1e-07);


tic
[kkappa,gap]=fminbnd(@(kkappa) ...
                    getKappaFunction(kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
                    nGridShocks,chi,upperBound,a,nAssets,...
                    ifLabor,nGridLabor,options),...
                    0.9704,0.99,options);
toc

% tic
% [~,laborParticipationRate]=getKappaFunction(kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
%     nGridShocks,chi,upperBound,a,nAssets,...
%     ifLabor,nGridLabor,options);
% toc                
[~,laborParticipationRate,kOverL]=getKappaFunction(kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor,options);

[~, laborParticipationRate,r,wage] = generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor)

table(gap,laborParticipationRate,r,wage,kkappa)

save('kkappaCalibration.mat')