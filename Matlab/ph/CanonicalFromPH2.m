%  [beta, B] = CanonicalFromPH2(alpha, A, prec)
%  
%  Returns the canonical form of an order-2 phase-type 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,2)
%    Initial vector of the phase-type distribution
%  A : matrix, shape (2,2)
%    Transient generator of the phase-type distribution
%  prec : double, optional
%    Numerical precision, default value is 1e-14
%  
%  Returns
%  -------
%  beta : matrix, shape (1,2)
%    The initial probability vector of the canonical form
%  B : matrix, shape (2,2)
%    Transient generator matrix of the canonical form
%  
%  Notes
%  -----
%  This procedure calculates 3 moments of the input and
%  calls 'PH2From3Moments'.

function [beta, B] = CanonicalFromPH2 (alpha, A, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha,A)
        error('CanonicalFromPH2: Input isn''t a valid ME distribution!');
    end

    if size(A,1)~=2
        error('CanonicalFromPH2: Dimension is not 2!');
    end

    [beta, B] = PH2From3Moments (MomentsFromME(alpha, A, 3), prec);
end

