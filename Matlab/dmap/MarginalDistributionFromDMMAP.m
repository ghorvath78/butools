%  [alpha, A] = MarginalDistributionFromDMMAP(D, precision)
%  
%  Returns the discrete phase type distributed marginal 
%  distribution of a discrete marked Markovian arrival 
%  process.
%  
%  Parameters
%  ----------
%  D : list/cell of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the DMMAP
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

function [alpha,A] = MarginalDistributionFromDMMAP(H)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMMAPRepresentation(H)
        error('MarginalDistributionFromDMMAP: Input isn''t a valid DMMAP representation');
    end

	[alpha,A] = MarginalDistributionFromDMRAP(H);
end