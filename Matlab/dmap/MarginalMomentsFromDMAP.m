%  moms = MarginalMomentsFromDMAP(D0, D1, K, precision)
%  
%  Returns the moments of the marginal distribution of a 
%  discrete Markovian arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the discrete Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the discrete Markovian arrival process
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

function moms = MarginalMomentsFromDMAP (D0, D1, K)

    if ~exist('K','var') || K==0
        K = 2*size(D0,1)-1;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMAPRepresentation(D0,D1)
        error('MarginalMomentsFromDMAP: input isn''t a valid DMAP representation!');
    end

    [alpha,A] = MarginalDistributionFromDMAP(D0,D1);
    moms = MomentsFromDPH(alpha,A,K);
end
