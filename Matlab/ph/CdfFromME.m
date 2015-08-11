%  cdf = CdfFromME(alpha, A, x)
%  
%  Returns the cummulative distribution function of a
%  matrix-exponential distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-exponential
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential
%      distribution.
%  x : vector of doubles
%      The cdf will be computed at these points
%  
%  Returns
%  -------
%  cdf : column vector of doubles
%      The values of the cdf at the corresponding "x" values

function cdf = CdfFromME (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('CdfFromME: Input isn''t a valid ME distribution!');
    end
    
    cdf = zeros(1,length(x));
    for i=1:length(x)
        cdf(i) = 1.0 - sum(alpha*expm(A*x(i)));
    end
end

