%  r = CheckProbMatrix(P, transient, prec)
%  
%  Checks if the matrix is a valid probability matrix: the 
%  matrix is a square matrix, the matrix has positive or 
%  zero off-diagonal elements, the rowsum of the matrix is 1.
%  
%  If "transient" is true, it checks if the matrix is a 
%  valid transient probability matrix: the matrix is a square
%  matrix, the matrix has positive or zero off-diagonal 
%  elements, the rowsum of the matrix is less than or equal
%  to 1, the maximum absolute eigenvalue is less than 1. 
%  
%  Parameters
%  ----------
%  P : matrix, shape (M,M)
%      The matrix to check.
%  transient : bool, optional
%      If true, the procedure checks if P is a transient 
%      probability matrix, otherwise it checks if it is
%      a valid probability matrix. The default value is 
%      false.
%  prec : double, optional
%      Entries with absolute value less than prec are 
%      considered to be zeros. The default value is 1e-14.
%      
%  Returns
%  -------
%  r : bool
%      The result of the check.

function r = CheckProbMatrix (P,transient,prec)

    global BuToolsVerbose;
    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end
   
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if ~exist('transient','var')
        transient = 0;
    end

    r = false;
    
    if size(P,1)~=size(P,2)
        if BuToolsVerbose
            fprintf ('CheckProbMatrix: the matrix is not a square matrix!\n');
        end
        return;
    end

    if min(min(P))<-prec
        if BuToolsVerbose
            fprintf ('CheckProbMatrix: the matrix has negative element (precision: %g)!\n', prec);
        end
        return;
    end

    if transient
        if any(sum(P,2)-1>size(P,2)*prec)
            if BuToolsVerbose
                fprintf ('CheckProbMatrix: The rowsum of the matrix (transient) is not less or equal than 1 (precision: %g)!\n', prec);
            end
            return;
        end

        if max(real(eig(P)))>=1-prec
            if BuToolsVerbose
                fprintf ('CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1 (precision: %g)!\n', prec);
            end
            return;
        end
    else
        if any(abs(sum(P,2)-1)>size(P,2)*prec)
            if BuToolsVerbose
                fprintf ('CheckProbMatrix: The rowsum of the matrix is not 1 (precision: %g)!\n', prec);
            end
            return;
        end
    end
    r = true;
end