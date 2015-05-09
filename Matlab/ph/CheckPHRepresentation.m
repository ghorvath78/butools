%  r = CheckPHRepresentation(alpha, A, prec)
%  
%  Checks if the given vector and matrix define a valid phase-
%  type representation.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      Initial vector of the phase-type distribution to check
%  A : matrix, shape (M,M)
%      Transient generator of the phase-type distribution to
%      check
%  prec : double, optional
%      Numerical precision. The default value is 1e-14.
%  
%  Returns
%  -------
%  r : bool
%      True, if vector alpha is a probability vector and matrix
%      A is a transient generator, and they have the same size.

function r = CheckPHRepresentation (alpha, A, prec)

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end

    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if length(alpha)~=size(A,1)
        if BuToolsVerbose
            fprintf ('CheckPHRepresentation:the vector and the matrix have different sizes!\n');
        end
        r = false;
        return;
    end

    if ~CheckGenerator(A,true,prec) || ~CheckProbVector(alpha,true,prec)
        r = false;
        return;
    end

    r = true;
end
