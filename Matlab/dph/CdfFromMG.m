%  cdf = CdfFromMG(alpha, A, x, prec)
%  
%  Returns the cummulative distribution function of a 
%  matrix-geometric distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-geometric distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-geometric 
%      distribution.
%  x : vector of non-negative integers
%      The density function will be computed at these points
%  prec : double, optional
%      Numerical precision to check if the input MG 
%      distribution is valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  cdf : column vector of doubles
%      The probabilities that the matrix-geometrically 
%      distributed random variable is less or equal to
%      the corresponding "x" values
%      

function cdf = CdfFromMG (alpha, A, x, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation(alpha, A, prec)
        error('CdfFromMG: Input isn''t a valid MG distribution!');
    end

    cdf = zeros(length(x),1);
    for i=1:length(x)
        cdf(i) = 1.0 - sum (alpha*(A^x(i)));
    end    
end

