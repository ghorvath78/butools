%  moms = MomentsFromME(alpha, A, K, prec)
%  
%  Returns the first K moments of a matrix-exponential
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-exponential
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential
%      distribution.
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments
%      are computed. The default value is K=0.
%  prec : double, optional
%      Numerical precision for checking the input.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  moms : row vector of doubles
%      The vector of moments.
%      

function moms = MomentsFromME (alpha, A, K)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('MomentsFromME: Input isn''t a valid ME representation!');
    end

    if ~exist('K','var') || K==0
        K = 2*length(alpha)-1;
    end

    moms = zeros(1,K);
    iA = inv(-A);
    for i=1:K
        moms(i) = factorial(i) * sum(alpha*iA^i);
    end
end
