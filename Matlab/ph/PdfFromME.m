%  pdf = PdfFromME(alpha, A, x, prec)
%  
%  Returns the probability density function of a matrix-
%  exponential distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-exponential
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential
%      distribution.
%  x : vector of doubles
%      The density function will be computed at these points
%  prec : double, optional
%      Numerical precision to check if the input ME 
%      distribution is valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  pdf : column vector of doubles
%      The values of the density function at the 
%      corresponding "x" values

function pdf = PdfFromME (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('PdfFromME: Input isn''t a valid ME distribution!');
    end

    pdf = zeros(1,length(x));
    for i=1:length(x)
       pdf(i) = sum(alpha*expm(A*x(i))*(-A));
    end
end
