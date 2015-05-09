%  r = CheckMMAPRepresentation(D, prec)
%  
%  Checks if the input matrixes define a continuous time MMAP.
%  
%  All matrices D0...DK must have the same size, D0 must be a 
%  transient generator matrix, D1 has only non-negative 
%  elements, and the rowsum of D0+D1+...+DK is 0 (up to the 
%  numerical precision).
%  
%  Parameters
%  ----------
%  D : list/cell of matrices, length(K)
%      The D0...DK matrices of the MMAP to check
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckMMAPRepresentation(H,prec)

    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    global BuToolsVerbose;

    if min(min(cell2mat(H(2:end)))) < -prec
        if BuToolsVerbose
            fprintf ('CheckMMAPRepresentation: Some of the matrices H1 ... HM have a negative element!\n');
        end
        r=false;
        return;
    end           

    r = CheckMAPRepresentation(H{1},SumMatrixList(H(2:end)),prec);
end
