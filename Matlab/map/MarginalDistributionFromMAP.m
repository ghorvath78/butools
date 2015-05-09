%  [alpha, A] = MarginalDistributionFromMAP(D0, D1, precision)
%  
%  Returns the phase type distributed marginal distribution
%  of a Markovian arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the Markovian arrival process
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the phase type 
%      distributed marginal distribution
%  A : matrix, shape (M,M)
%      The transient generator of the phase type distributed
%      marginal distribution    

function [alpha,A] = MarginalDistributionFromMAP (D0, D1)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1)
        error('MarginalDistributionFromMAP: input isn''t a valid MAP representation!');
    end

    [alpha,A] = MarginalDistributionFromRAP(D0,D1);
end
