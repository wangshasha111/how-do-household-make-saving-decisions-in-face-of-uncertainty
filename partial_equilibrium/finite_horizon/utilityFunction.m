function [utility] = utilityFunction(consumption, ggamma)
% CRRA utility function
%     if consumption <= 0
%         utility = - 1000000000;
%     else
        if ggamma == 1
            utility = log(consumption);
        else
            utility = (consumption.^(1-ggama)-1)/(1-ggama);
        end
%     end
end