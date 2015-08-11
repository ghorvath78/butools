%  r = CheckGenerator(Q, transient, prec)
%  
%  Checks if the matrix is a valid generator matrix: the 
%  matrix is a square matrix, the matrix has positive or 
%  zero off-diagonal elements, the diagonal of the matrix 
%  is negative, the rowsum of the matrix is 0.
%  
%  If the "transient" parameter is set to false, it checks 
%  if the real part of the maximum absolute eigenvalue is 
%  less than zero and the rowsum is equal or less than 0. 
%  
%  Parameters
%  ----------
%  Q : matrix, shape (M,M)
%      The generator to check.
%  transient : bool, optional
%      If true, the procedure checks if Q is a transient 
%      generator, otherwise it checks if it is a valid 
%      generator. The default value is false.
%  prec : double, optional
%      Entries with absolute value less than prec are 
%      considered to be zeros. The default value is 1e-14.
%      
%  Returns
%  -------
%  r : bool
%      The result of the check.

function r = CheckGenerator (Q,transient,prec)

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
        transient = false;
    end

    r = false;

    if size(Q,1)~=size(Q,2)
        if BuToolsVerbose
            fprintf ('CheckGenerator: Generator is not a square matrix!\n');
        end
        return;
    end

    if any(diag(Q)>=prec)
        if BuToolsVerbose
            fprintf ('CheckGenerator: The diagonal of the generator is not negative (precision: %g)!\n', prec);
        end
        return;
    end

    N = size(Q,1);
    odQ = Q<-prec;
    for i=1:N
        odQ(i,i) = 0;
    end

    if sum(any(odQ))>0
        if BuToolsVerbose
            fprintf ('CheckGenerator: The generator has negative off-diagonal element (precision: %g)!\n', prec);
        end
        return;
    end

    if transient
        if max(sum(Q,2))>prec
            if BuToolsVerbose
                fprintf ('CheckGenerator: The rowsum of the transient generator is greater than 0 (precision: %g)!\n', prec);
            end
            return;
        end

        if max(real(eig(Q)))>=prec
            if BuToolsVerbose
                fprintf ('CheckGenerator: The transient generator has non-negative eigenvalue (precision: %g)!\n', prec);
            end
            return;
        end
    else
        if any(abs(sum(Q,2))>prec)
            if BuToolsVerbose
                fprintf ('CheckGenerator: The rowsum of the generator is not 0 (precision: %g)!\n', prec);
            end
            return;
        end
    end
    r = true;
end