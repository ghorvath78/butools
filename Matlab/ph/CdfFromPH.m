%  cdf = CdfFromPH(alpha, A, x)
%  
%  Returns the cummulative distribution function of a
%  continuous phase-type distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the phase-type
%      distribution.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the phase-type
%      distribution.
%  x : vector of doubles
%      The cdf will be computed at these points
%  
%  Returns
%  -------
%  cdf : column vector of doubles
%      The values of the cdf at the corresponding "x" values

function cdf = CdfFromPH (alpha, A, x)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckPHRepresentation(alpha, A)
        error('CdfFromPH: Input isn''t a valid PH distribution!');
    end

    cdf = CdfFromME (alpha, A, x);
end

