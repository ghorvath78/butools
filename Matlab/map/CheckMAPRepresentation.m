%  r = CheckMAPRepresentation(D0, D1, prec)
%  
%  Checks if the input matrixes define a continuous time MAP.
%  
%  Matrices D0 and D1 must have the same size, D0 must be a 
%  transient generator matrix, D1 has only non-negative 
%  elements, and the rowsum of D0+D1 is 0 (up to the numerical
%  precision).
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the MAP to check
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the MAP to check
%  prec : double, optional
%      Numerical precision, the default value is 1e-14
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckMAPRepresentation (D0, D1, prec)

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if ~CheckGenerator(D0,1,prec)
        r = false;
        return;
    end

    if size(D0)~=size(D1)
        if BuToolsVerbose
            fprintf ('CheckMAPRepresentation: D0 and D1 have different sizes!\n');
        end
        r = false;
        return;
    end

    if min(min(D1))<-prec
        if BuToolsVerbose
            fprintf ('CheckMAPRepresentation: D1 has negative element (precision: %g)!\n', prec);
        end
        r = false;
        return;
    end

    if any(abs(sum(D0+D1,2))>prec)
        if BuToolsVerbose
            fprintf ('CheckMAPRepresentation: The rowsum of D0+D1 is not 0 (precision: %g)!\n', prec);
        end
        r = false;
        return;
    end

    r = true;
end
