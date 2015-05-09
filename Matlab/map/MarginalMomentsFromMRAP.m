%  moms = MarginalMomentsFromMRAP(H, K, precision)
%  
%  Returns the moments of the marginal distribution of a 
%  marked rational arrival process.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the MRAP
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

function moms = MarginalMomentsFromMRAP (H,K)

    if ~exist('K','var') || K==0
        K = 2*size(H{1},1)-1;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   
    
    if BuToolsCheckInput && ~CheckMRAPRepresentation (H)
        error('MarginalMomentsFromMRAP: Input isn''t a valid MRAP representation!');
    end

    [alpha,A] = MarginalDistributionFromMRAP(H);
    moms = MomentsFromME(alpha,A,K);
end
