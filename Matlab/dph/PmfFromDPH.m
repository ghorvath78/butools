%  pmf = PmfFromDPH(alpha, A, x, prec)
%  
%  Returns the probability mass function of a discrete
%  phase-type distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial probability vector of the discrete phase-
%      type distribution. The sum of the entries of pi0 is
%      less or equal to 1.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the discrete phase-
%      type distribution.
%  x : vector of non-negative integers
%      The density function will be computed at these points
%  prec : double, optional
%      Numerical precision to check if the input DPH 
%      distribution is valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  pmf : column vector of doubles
%      The probabilities that the discrete phase type
%      distributed random variable takes the corresponding
%      "x" values
%      

function pmf = PmfFromDPH (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDPHRepresentation(alpha, A)
        error('PmfFromDPH: Input isn''t a valid DPH distribution!');
    end

    pmf = PmfFromMG (alpha, A, x);
end

