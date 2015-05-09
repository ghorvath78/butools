%  pdf = PdfFromPH(alpha, A, x, prec)
%  
%  Returns the probability density function of a continuous
%  phase-type distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial probability vector of the phase-type
%      distribution.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the phase-type
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

function pdf = PdfFromPH (alpha, A, x)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckPHRepresentation(alpha, A)
        error('PdfFromPH: Input isn''t a valid PH distribution!');
    end

    pdf = PdfFromME (alpha, A, x);
end
