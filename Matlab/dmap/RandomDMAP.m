%  [D0, D1] = RandomDMAP(order, mean, zeroEntries, maxTrials, prec)
%  
%  Returns a random disctere Markovian arrival process.
%  
%  Parameters
%  ----------
%  order : int
%      The size of the DMAP
%  mean : double, optional
%      The mean inter-arrival times of the DMAP
%  zeroEntries : int, optional
%      The number of zero entries in the D0 and D1 matrices
%  maxTrials : int, optional
%      The maximum number of trials to find a proper DMAP 
%      (that has an irreducible phase process and none of 
%      its parameters is all-zero)
%  prec : double, optional
%      Numerical precision for checking the irreducibility.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  D0 : vector, shape (1,M)
%      The D0 matrix of the DMAP
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the DMAP
%  
%  Notes
%  -----
%  If it fails, try to increase the 'maxTrials' parameter,
%  or/and the 'mean' parameter.

function [d0,d1] = RandomDMAP(order, mean, zeroEntries, maxTrials, prec)
% RandomDMAP [ order, zeroEntries[0], prec[10^-14], maxTrials[100] ]
%  -> [ matrix0, matrix1 ] : 
%     Generates a random discrete time Markovian arrival process of the 
%     given order. The obtained representation containes 'zeroEntries' zeros.
%     If it fails after 'maxTrials' trials, then it decreases the number of
%     zero entries. It prints a message, if the found representation contains
%     less zeros than given. 'prec' is the numerical precision.

    if ~exist('zeroEntries','var')
        zeroEntries = 0;
    end

    if ~exist('prec','var')
        prec = 1e-7;
    end

    if ~exist('mean','var')
        mean = 10;
    end
    
    if ~exist('maxTrials','var')
        maxTrials = 1000;
    end

    D = RandomDMMAP(order, 1, mean, zeroEntries, maxTrials, prec);
    d0=D{1};
    d1=D{2};
end