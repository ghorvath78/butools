%  r = CheckMEPositiveDensity(alpha, A, maxSize, prec)
%  
%  Checks if the given matrix-exponential distribution has 
%  positive density.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      Initial vector of the matrix-exponential distribution 
%      to check
%  A : matrix, shape (M,M)
%      Matrix parameter of the matrix-exponential distribution
%      to check
%  maxSize : int, optional
%      The procedure tries to transform the ME distribution
%      to phase-type up to order maxSize. The default value
%      is 100.
%  prec : double, optional
%      Numerical precision. The default value is 1e-14.
%  
%  Returns
%  -------
%  r : bool
%      True, if the given matrix-exponential distribution has
%      a positive density
%  
%  Notes
%  -----
%  This procedure calls MonocyclicPHFromME, and can be time 
%  consuming. 

function r = CheckMEPositiveDensity (alpha, A, maxSize, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('maxSize','var')
        maxSize = 100;
    end

    try
        [beta, B] = MonocyclicPHFromME (alpha, A, maxSize, prec);
        r = CheckMERepresentation (alpha, A, prec);
    catch err
        r = false;
    end
end

