function [utility] = utilityFunction(consumption, ggamma)
% CRRA utility function
    if consumption <= 0
        utility = - inf;
    else
        if ggamma == 1
            utility = log(consumption);
        else
            utility = ((consumption.^(1-ggamma))-1)/(1-ggamma);
        end
    end
return