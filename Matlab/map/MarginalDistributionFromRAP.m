%  [alpha, A] = MarginalDistributionFromRAP(H0, H1, precision)
%  
%  Returns the phase type distributed marginal distribution
%  of a rational arrival process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
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

function [alpha,A] = MarginalDistributionFromRAP (H0, H1)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckRAPRepresentation(H0,H1)
        error('MarginalDistributionFromRAP: Input isn''t a valid RAP representation!');
    end

    alpha = DRPSolve(inv(-H0)*H1);
    A = H0;
end

