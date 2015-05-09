%  moms = MomentsFromPH(alpha, A, K, prec)
%  
%  Returns the first K moments of a continuous phase-type
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial probability vector of the phase-type
%      distribution.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the phase-type
%      distribution.
%  K : int, optional
%      Number of moments to compute. If K=0, 2*M-1 moments
%      are computed. The default value is K=0.
%  prec : double, optional
%      Numerical precision for checking the input.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  moms : row vector of doubles
%      The vector of moments.

function moms = MomentsFromPH (alpha, A, K)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckPHRepresentation(alpha, A)
        error('MomentsFromPH: Input isn''t a valid PH representation!');
    end

    if ~exist('K','var')
        K = 0;
    end

    moms = MomentsFromME (alpha, A, K);
end
