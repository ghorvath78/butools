%  r = CheckDMMAPRepresentation(D, prec)
%  
%  Checks if the input matrixes define a discrete time MMAP.
%  
%  All matrices D0...DK must have the same size, D0 must be a 
%  transient probability matrix, D1 has only non-negative 
%  elements, and the rowsum of D0+D1+...+DK is 1 (up to the 
%  numerical precision).
%  
%  Parameters
%  ----------
%  D : list/cell of matrices, length(K)
%      The D0...DK matrices of the DMMAP to check
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckDMMAPRepresentation(D,prec)
% CheckDMMAPRepresentation [ (vector of D0, D1 .. DM), prec[10^-14] ] : 
%     Checks if the input matrixes define a discrete time MMAP: D0 and
% 	  the sum of D1..DM define a DMAP. 'prec' is the numerical precision.

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if min(min(cell2mat(D(2:end)))) < -prec
        if BuToolsVerbose
            fprintf ('CheckDMMAPRepresentation: One of the matrices D1 ... DM has a negative element!\n');
        end
        r=0;
        return;
    end

    r = CheckDMAPRepresentation(D{1},SumMatrixList(D(2:end)),prec);

end