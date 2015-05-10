%  [alpha, A] = MarginalDistributionFromDRAP(H0, H1, precision)
%  
%  Returns the matrix geometrically distributed marginal 
%  distribution of a discrete rational arrival process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the discrete rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the discrete rational arrival process
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

function [alpha,A] = MarginalDistributionFromDRAP( H0, H1)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDRAPRepresentation(H0,H1)
        error('MarginalDistributionFromDRAP: Input isn''t a valid DRAP representation!');
    end

    alpha = DRPSolve(inv(eye(size(H0,1))-H0)*H1);
    A=H0;
end