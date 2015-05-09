%  [beta, B] = CanonicalFromPH3(alpha, A, prec)
%  
%  Returns the canonical form of an order-3 phase-type 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,3)
%      Initial vector of the phase-type distribution
%  A : matrix, shape (3,3)
%      Transient generator of the phase-type distribution
%  prec : double, optional
%    Numerical precision, default value is 1e-14
%  
%  Returns
%  -------
%  beta : matrix, shape (1,3)
%    The initial probability vector of the canonical form
%  B : matrix, shape (3,3)
%    Transient generator matrix of the canonical form
%  
%  Notes
%  -----
%  This procedure calculates 5 moments of the input and
%  calls 'PH3From5Moments'.

function [beta, B] = CanonicalFromPH3 (alpha, A, prec)

    if ~exist('prec','var')
        prec = 1e-10;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha,A)
        error('CanonicalFromPH3: Input isn''t a valid ME distribution!');
    end

    if size(A,1)~=3
        error('CanonicalFromPH3: Dimension is not 3!');
    end

    [beta, B] = PH3From5Moments (MomentsFromME(alpha, A, 5), prec);
end

