function [laborParticipationRateGap,laborParticipationRate,kOverL] = getKappaFunction(kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
                    nGridShocks,chi,upperBound,a,nAssets,...
                    ifLabor,nGridLabor,options)
                
[kOverL,diff]=fminbnd(@(kOverL) ...
                    generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
                    nGridShocks,chi,upperBound,a,nAssets,...
                    ifLabor,nGridLabor),...
                    0.000001,10,options);
                
[~,laborParticipationRate] = generalEqEndoLaborFunction(kOverL,kkappa, ggamma, ddelta, ssigmaY,bbeta,aalphaK,depreciation,TFP,...
    nGridShocks,chi,upperBound,a,nAssets,...
    ifLabor,nGridLabor)

    table(kOverL,diff,laborParticipationRate)
    laborParticipationRateGap = abs(laborParticipationRate-0.8);
    return
%                 vKOverL(iKappa) = kOverL;
%                 vDiff(iKappa) = diff;
%                 vLaborParticipationRate(iKappa) = laborParticipationRate;
                
% end
% 
% figure
% plot(vKappa,vKOverL)
% 
% figure
% plot(vKappa,vDiff)
% 
% figure
% plot(vKappa,vLaborParticipationRate)
% 
% [v,ind] = min(vLaborParticipationRate-0.8);
% 
% kkappaCalibrated = vKappa(ind);
% 
% save('kkappaCalibration')