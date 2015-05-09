%  [alpha, A] = MarginalDistributionFromMRAP(H, precision)
%  
%  Returns the phase type distributed marginal distribution
%  of a marked rational arrival process.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the MRAP
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix exponentially
%      distributed marginal
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix exponentially
%      distributed marginal    

function [alpha,A] = MarginalDistributionFromMRAP (H)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   
    
    if BuToolsCheckInput && ~CheckMRAPRepresentation (H)
        error('MarginalDistributionFromMRAP: Input isn''t a valid MRAP representation!');
    end

    alpha = DRPSolve(inv(-H{1})*SumMatrixList(H(2:end)));
    A=H{1};
end
