%% Question 3.2 
% Reproduce table II from Aiyagari (1994)

% households
rrho = 0.04; % discount rate
bbeta = 1/(1+rrho); % discount factor

% firms
aalphaK = 0.36;            % capital share
depreciation = 0.08;     % depreciation rate
labourEq = 1; % Normalize labor to 1
TFP = 1; % total production factor

% income shocks
nGridShocks = 31;
chi = 0; % tightness of borrowing constraint, between 0 (no borrowing at all) and 1 ( just to avoid Ponzi scheme)
upperBound = 100 ; % TRY DIFFERENT VALUES % ideally 2 times the steady state?

% Grid of assets
a = 0.2;    % spacing of the grid (exponential)
nAssets = 1000;

% Different Parameters
vDelta = [0,0.3,0.6,0.9]; % income process persistence
vSigmaY = [0.2,0.4];             % income process variance
vGamma = [1,3,5];        % 1/EIS

% Optimization
options = optimset('Display','Iter','TolX',1e-07);  
mEquilibriumRate = zeros(length(vDelta),length(vGamma),length(vSigmaY));
mDifference = zeros(length(vDelta),length(vGamma),length(vSigmaY));

% take initial guesses from Aiyagari (1994) table II
mInitialGuess = [4.15, 4.1365, 4.0912, 3.9305;
                 4.1456, 4.0432, 3.8767, 3.2903;
                 4.0858, 3.9054, 3.5857, 2.4260];
% mInitialGuess = [4.0649, 3.9554, 3.7567, 3.3054;
%                  3.7816, 3.4188, 2.7835, 1.2894;
%                  3.4177, 2.8032, 1.8070, -0.3456];
             
mInitialGuess = mInitialGuess/100;

iteration = 0;
totalIteration = length(vDelta) * length(vSigmaY) * length(vGamma);
          
% tic
for iSigmaY = 1 : length(vSigmaY)
    
    ssigmaY = vSigmaY(iSigmaY);
    for iDelta=1:length(vDelta)
        ddelta = vDelta(iDelta);
        
        for iGamma=1:length(vGamma)
            ggamma = vGamma(iGamma);

    %         parameters.delta = vDelta(iDelta);
    %         parameters.gamma = vGamma(iGamma);
            guess = mInitialGuess(iGamma,iDelta);

            myfn = @(r) generalEqFunction(r,ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,labourEq,TFP,...
                nGridShocks,chi,upperBound,a,nAssets);

            tic                             
            [mEquilibriumRate(iDelta,iGamma,iSigmaY),mDifference(iDelta,iGamma,iSigmaY)] = fzero(myfn,guess,options);
            table([mEquilibriumRate(iDelta,iGamma,iSigmaY),mDifference(iDelta,iGamma,iSigmaY),ssigmaY,ggamma,ddelta])
            toc
            
            iteration = iteration +1;
            fprintf('Number of Iteration: %2.0f\n Progress of Iteration: %2.0f%%\n',iteration,round(iteration/totalIteration,2)*100);

        end
        iteration = iteration +1;
        fprintf('Number of Iteration: %2.0f\n Progress of Iteration: %2.0f%%\n',iteration,round(iteration/totalIteration,2)*100);        
        
    end
    
    iteration = iteration +1;
    fprintf('Number of Iteration: %2.0f\n Progress of Iteration: %2.0f%%\n',iteration,round(iteration/totalIteration,2)*100);
end
% toc

save('AiyagariTable')

% I tried different tolerance level of the fzero function, but found out
% that the accuracy of the calibration of r, i.e., how small the gap between
% supply and demand can be, depends HEAVILY on how fine the asset grid is,
% NOT on the range of the grid, NOT on how fine the income shock grid is.
% Once you have many grid points for assets, make sure that you set the
% tolerance level to 1e-7 so that the resulted accuracy could be around
% 1e-4.
 