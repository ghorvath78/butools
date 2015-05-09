%  moms = MarginalMomentsFromRAP(H0, H1, K, precision)
%  
%  Returns the moments of the marginal distribution of a 
%  rational arrival process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments
%      are computed. The default value is K=0.
%  precision : double, optional
%      Numerical precision for checking if the input is valid.
%      The default value is 1e-14
%  
%  Returns
%  -------
%  moms : row vector of doubles, length K
%      The vector of moments.

function moms = MarginalMomentsFromRAP (H0, H1, K)

    if ~exist('K','var') || K==0
        K = 2*size(H0,1)-1;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckRAPRepresentation(H0,H1)
        error('MarginalMomentsFromRAP: Input isn''t a valid RAP representation!');
    end

    [alpha,A] = MarginalDistributionFromRAP(H0,H1);
    moms = MomentsFromME(alpha,A,K);
end

