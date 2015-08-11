%  [x,y] = IntervalPdfFromPH(alpha, A, intBounds, prec)
%  
%  Returns the approximate probability density function of a
%  continuous phase-type distribution, based on the 
%  probability of falling into intervals.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial probability vector of the phase-type
%      distribution.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the phase-type
%      distribution.
%  intBounds : vector, shape (K)
%      The array of interval boundaries. The pdf is the
%      probability of falling into an interval divided by
%      the interval length. 
%      If the size of intBounds is K, the size of the result is K-1.
%  prec : double, optional
%      Numerical precision to check if the input is a valid
%      phase-type distribution. The default value is 1e-14
%  
%  Returns
%  -------
%  x : matrix of doubles, shape(K-1,1)
%      The points at which the pdf is computed. It holds the center of the 
%      intervals defined by intBounds.
%  y : matrix of doubles, shape(K-1,1)
%      The values of the density function at the corresponding "x" values
%  
%  Notes
%  -----
%  This method is more suitable for comparisons with empirical
%  density functions than the exact one (given by PdfFromPH).

function [x, y] = IntervalPdfFromPH (alpha, A, intBounds)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckPHRepresentation(alpha, A)
        error('IntervalPdfFromPH: Input isn''t a valid PH distribution!');
    end
    
    K = length(intBounds);
    x = reshape ((intBounds(2:end) + intBounds(1:end-1)) / 2, K-1,1);
    y = zeros(1,K-1);
    for i=1:K-1
        y(i) = (sum(alpha*expm(A*intBounds(i))) - sum(alpha*expm(A*intBounds(i+1))))/(intBounds(i+1)-intBounds(i));
    end   
end
