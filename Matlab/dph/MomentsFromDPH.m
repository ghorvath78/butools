%  moms = MomentsFromDPH(alpha, A, K, prec)
%  
%  Returns the first K moments of a discrete phase-type
%  distribution.
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
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments
%      are computed. The default value is 0.
%  prec : double, optional
%      Numerical precision for checking the input.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  moms : row vector of doubles
%      The vector of moments.
%      

function moms = MomentsFromDPH (alpha, A, K)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDPHRepresentation(alpha, A)
        error('MomentsFromDPH: Input isn''t a valid PH representation!');
    end

    if ~exist('K','var')
        K = 0;
    end
    
    moms = MomentsFromMG (alpha, A, K);
end

