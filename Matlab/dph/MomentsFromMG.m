%  moms = MomentsFromMG(alpha, A, K, prec)
%  
%  Returns the first K moments of a matrix geometric 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-geometric distribution.
%      The sum of the entries of alpha is less or equal to 1.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-geometric 
%      distribution.
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments are
%      computed. The default value is 0.
%  prec : double, optional
%      Numerical precision for checking the input.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  moms : row vector of doubles
%      The vector of moments.
%      

function moms = MomentsFromMG (alpha, A, K)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation (alpha, A)
        error('MomentsFromMG: Input isn''t a valid MG representation!');
    end

    if ~exist('K','var') || K==0
        K = 2*length(alpha) - 1;
    end

    fmoms = zeros(1,K);
    iA = inv(eye(size(A))-A);
    for i=1:K
        fmoms(i) = factorial(i) * sum(alpha*iA^i*A^(i-1));
    end
    moms = MomsFromFactorialMoms (fmoms);
end

