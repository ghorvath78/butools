%  pi = FluidStationaryDistr (mass0, ini, K, clo, x)
%  
%  Returns the stationary distribution of a Markovian 
%  fluid model at the given points.
%  
%  Parameters
%  ----------
%  mass0 : matrix, shape (1,Np+Nm)
%      The stationary probability vector of zero level
%  ini : matrix, shape (1,Np)
%      The initial vector of the stationary density
%  K : matrix, shape (Np,Np)
%      The matrix parameter of the stationary density
%  clo : matrix, shape (Np,Np+Nm)
%      The closing matrix of the stationary density
%  x : vector, length (K)
%      The distribution function is computed at these 
%      points.
%  
%  Returns
%  -------
%  pi : matrix, shape (K,Nm+Np)
%      The ith row of pi is the probability that the fluid
%      level is less than or equal to x(i), while being in
%      different states of the background process.

function y = FluidStationaryDistr (mass0, ini, K, clo, x)

    m = size(clo,2);
    y = zeros(length(x),m);
    closing = inv(-K)*clo;
    for i=1:length(x)
        y(i,:) = mass0 + ini*(eye(size(K))-expm(K*x(i)))*closing;
    end
end

