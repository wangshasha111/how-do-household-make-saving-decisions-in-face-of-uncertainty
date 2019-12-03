function [Z,Zprob] = tauchen(N,mu,rho,sigma,m)

% Purpose:    Finds a N-Markov chain that approximates the AR(1) process
%
%                z(t+1) = (1-rho)*mu+rho*z(t)+sigma*eps(t+1)
%
% Format:     [Z,Zprob] = tauchen(N,mu,rho,sigma,m)
%
% Input:      N       scalar, number of nodes for Z
%             mu      scalar, unconditional mean of process
%             rho     scalar
%             sigma   scalar, std. dev. of epsilons
%             m       max +- std. devs, 3 is suggested.
%
% Output:     Z       N*1 vector, nodes for Z
%             Zprob   N*N matrix, transition probabilities
%
%    Martin Floden, Fall 1996
%    Modified by Dario Caldara and Jesus Fernandez-Villaverde, Fall 2007
%     Modified by Rodrigo Morales 2019 (added case N = 1);
%       case N = 1 is important for quick checks and stationary case.
%
%    This procedure implements George Tauchen's algorithm described in Ec. Letters 20 (1986) 177-181.
if N > 1
    Z     = zeros(N,1);
    Zprob = zeros(N,N);
    a     = (1-rho)*mu;

    Z(N)  = m*sqrt(sigma^2/(1-rho^2));
    Z(1)  = -Z(N);
    zstep = (Z(N)-Z(1))/(N-1);

    for i=2:(N-1)
        Z(i) = Z(1)+zstep*(i-1);
    end 

    Z = Z+a/(1-rho);

    for j = 1:N
        for k = 1:N
            if k == 1
                Zprob(j,k) = normcdf((Z(1)-a-rho*Z(j)+zstep/2)/sigma);
            elseif k == N
                Zprob(j,k) = 1-normcdf((Z(N)-a-rho*Z(j)-zstep/2)/ sigma);
            else
                Zprob(j,k) = normcdf((Z(k)-a-rho*Z(j)+zstep/2)/sigma)-...
                             normcdf((Z(k)-a-rho*Z(j)-zstep/2)/sigma);
            end
        end
    end

    Z=Z';
else
    Z = mu;
    Zprob = 1;
end