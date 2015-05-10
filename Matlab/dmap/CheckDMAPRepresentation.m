%  r = CheckDMAPRepresentation(D0, D1, prec)
%  
%  Checks if the input matrixes define a discrete time MAP.
%  
%  Matrices D0 and D1 must have the same size, D0 must be a 
%  transient probability matrix, D1 has only non-negative
%  elements, and the rowsum of D0+D1 is 1 (up to the numerical
%  precision).
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the DMAP to check
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the DMAP to check
%  prec : double, optional
%      Numerical precision, the default value is 1e-14
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r=CheckDMAPRepresentation(d0, d1, prec)

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if ~CheckProbMatrix(d0,1,prec)
        if BuToolsVerbose
            fprintf('CheckDMAPRepresentation: D0 isn''t a transient probability matrix!\n');
        end
        r=false;
        return;
    end

    if size(d0,1) ~= size(d1,1) || size(d0,2) ~= size(d1,2)
        if BuToolsVerbose
            fprintf('CheckDMAPRepresentation: D0 and D1 have different sizes!\n');
        end
        r=false;
        return;
    end

    if min(min([d0,d1])) < -prec
        if BuToolsVerbose
            fprintf('CheckDMAPRepresentation: One of the matrices has negative element!\n');
        end
        r=false;
        return;
    end

    if any(abs(sum(d0+d1,2)-1) > prec)
        if BuToolsVerbose
            fprintf('CheckDMAPRepresentation: A rowsum of matrix0+matrix1 is not 1 (at precision %g)!\n',prec);
        end
        r=false;
        return;
    end

    r=true;

end