%  pmf = PmfFromMG(alpha, A, x, prec)
%  
%  Returns the probability mass function of a matrix-
%  geometric distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-geometric
%      distribution. The sum of the entries of pi0 is less
%      or equal to 1.
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
%  pmf : column vector of doubles
%      The probabilities that the matrix-geometrically 
%      distributed random variable takes the corresponding "x"
%      values
%      

function pmf = PmfFromMG (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation(alpha, A)
        error('PmfFromMG: Input isn''t a valid MG distribution!');
    end

    a = 1.0 - sum(A,2);
    pmf = zeros(1,length(x));
    for i=1:length(x)
        if x(i)==0
            pmf(i) = 1.0 - sum(alpha);
        else
            pmf(i) = sum(alpha*(A^(x(i)-1))*a);
        end
    end
end
