%% Research Project 
% Economics 712
% Fall 2019
% Professor Dirk Krueger

% Value function iteration partially taken from "Basic RBC model with full depreciation" by Jesus Fernandez-Villaverde Haverford, July 3, 2013

% 2019-11-6 22:06:22
% 2019-12-3

clear;
close all;

% Part I Partial Equilibrium

%% model and parameters
ggamma = 1; % inverse of intertemporal elasticity of substitution
ddelta = 0.7; % persistence of the income process
ssigmaY = [0.2,0.25]; % income process variance
ssigmaError=ssigmaY*(sqrt(1-ddelta^2)); % variance of the error term of the income process
ssigmaErrorLow = min(ssigmaError);
ssigmaErrorHigh = max(ssigmaError);

rrho = 0.1; % discount rate
bbeta = 1/(1+rrho); % discount factor
r = 0.02;% interest rate

% income shocks
nGridShocks = 51;
% **************************************************** TRY DIFFERENT VALUES ****************************************************
[vIncomeShocks, mTransition]=rouwenhorstFunction(ddelta,ssigmaErrorLow,nGridShocks); % TRY DIFFERENT VALUES of ssigmaError, low and high
vIncomeShocks = exp(vIncomeShocks)';

% assets
AssetLimit = min(vIncomeShocks); % worst case scenario of income process
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
minK = -chi*AssetLimit/r; % natural borrowing limit

% **************************************************** TRY DIFFERENT VALUES ****************************************************
upperBound = 4 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?
maxK = upperBound; 
a = 0.2;    % spacing of the grid (exponential)
% **************************************************** TRY DIFFERENT VALUES ****************************************************
nAssets = 500;
if chi == 0
    temp = exp(linspace(log(minK+0.000001)*a,log(maxK)*a,nAssets))';
    temp_01 = (temp-min(temp))/(max(temp)-min(temp)); % re-scale to be between 0 and 1
    vGridAsset = temp_01*(maxK-minK) + minK; % re-scale to be between limits
    clear temp temp_01;
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset));% life starting point of asset
    vGridAsset(iAsset0)=0;
else
    vGridAsset = curvspaceFunction(minK, maxK, nAssets ,2)';
    nGridAsset = length(vGridAsset);
    [~, iAsset0]=min(abs(vGridAsset)); % life starting point of asset
end

%% exercise 7 Finite Horizon with deterministic retirement benefit
nSimulations = 50000; % number of simulated paths M
nPeriods = 61; % horizon

% import data
% readmatrix