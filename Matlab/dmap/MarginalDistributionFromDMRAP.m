%  [alpha, A] = MarginalDistributionFromDMRAP(H, precision)
%  
%  Returns the matrix geometrically distributed marginal 
%  distribution of a discrete marked rational arrival 
%  process.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the DMRAP
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix geometrically
%      distributed marginal distribution
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix geometrically
%      distributed marginal distribution    

function [alpha,A] = MarginalDistributionFromDMRAP(H)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMRAPRepresentation(H)
        error('MarginalDistributionFromDMRAP: Input isn''t a valid DMRAP representation');
    end

    alpha=DRPSolve(inv(eye(size(H{1},1))-H{1})*SumMatrixList(H(2:end)));
    A=H{1};

end