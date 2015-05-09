%  [D0, D1] = RandomMAP(order, mean, zeroEntries, maxTrials, prec)
%  
%  Returns a random Markovian arrival process with given mean 
%  value.
%  
%  Parameters
%  ----------
%  order : int
%      The size of the MAP
%  mean : double, optional
%      The mean inter-arrival times of the MAP
%  zeroEntries : int, optional
%      The number of zero entries in the D0 and D1 matrices
%  maxTrials : int, optional
%      The maximum number of trials to find a proper MAP 
%      (that has an irreducible phase process and none of 
%      its parameters is all-zero)
%  prec : double, optional
%      Numerical precision for checking the irreducibility.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  D0 : vector, shape (1,M)
%      The D0 matrix of the MAP
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the MAP

function [D0,D1] = RandomMAP(order, mean, zeroEntries, maxTrials, prec)

    if ~exist('zeroEntries','var')
        zeroEntries = 0;
    end

    if ~exist('mean','var')
        mean = 1;
    end

    if ~exist('prec','var')
        prec = 1e-7;
    end

    if ~exist('maxTrials','var')
        maxTrials = 1000;
    end

    D = RandomMMAP(order, 1, mean, zeroEntries, maxTrials, prec);
    D0=D{1};
    D1=D{2};
end
