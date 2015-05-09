%  cdf = CdfFromDPH(alpha, A, x, prec)
%  
%  Returns the cummulative distribution function of a 
%  discrete phase-type distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the discrete phase-
%      type distribution.
%  A : matrix, shape (M,M)
%      The transition probability  matrix of the discrete phase-
%      type distribution.
%  x : vector of non-negative integers
%      The cdf will be computed at these points
%  prec : double, optional
%      Numerical precision to check if the input DPH 
%      distribution is valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  cdf : column vector of doubles
%      The probabilities that the discrete phase type 
%      distributed random variable is less or equal to the
%      corresponding "x" values
%      

function cdf = CdfFromDPH (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDPHRepresentation(alpha, A)
        error('CdfFromDPH: Input isn''t a valid PH distribution!');
    end
    
    cdf = CdfFromMG (alpha, A, x);
end

