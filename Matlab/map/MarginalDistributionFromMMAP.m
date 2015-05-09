%  [alpha, A] = MarginalDistributionFromMMAP(D, precision)
%  
%  Returns the phase type distributed marginal distribution
%  of a marked Markovian arrival process.
%  
%  Parameters
%  ----------
%  D : list/cell of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the MMAP
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

function [alpha,A] = MarginalDistributionFromMMAP (D)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMMAPRepresentation(D)
        error('MarginalDistributionFromMMAP: Input isn''t a valid MMAP representation');
    end

    [alpha,A] = MarginalDistributionFromMRAP(D);
end
