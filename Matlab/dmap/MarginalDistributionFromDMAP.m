%  [alpha, A] = MarginalDistributionFromDMAP(D0, D1, precision)
%  
%  Returns the discrete phase type distributed marginal 
%  distribution of a discrete Markovian arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the discrete Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the discrete Markovian arrival process
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the discrete phase 
%      type distributed marginal distribution
%  A : matrix, shape (M,M)
%      The transient generator of the discrete phase type 
%      distributed marginal distribution    

function [alpha,A] = MarginalDistributionFromDMAP( D0, D1)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMAPRepresentation(D0,D1)
        error('MarginalDistributionFromDMAP: input isn''t a valid DMAP representation!');
    end

    [alpha,A] = MarginalDistributionFromDRAP(D0,D1);
end