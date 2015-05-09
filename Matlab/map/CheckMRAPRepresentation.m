%  r = CheckMRAPRepresentation(H, prec)
%  
%  Checks if the input matrixes define a continuous time MRAP.
%  
%  All matrices H0...HK must have the same size, the dominant
%  eigenvalue of H0 is negative and real, and the rowsum of 
%  H0+H1+...+HK is 0 (up to the numerical precision).
%  
%  Parameters
%  ----------
%  H : list/cell of matrices, length(K)
%      The H0...HK matrices of the MRAP to check
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckMRAPRepresentation(H,prec)

    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    r = CheckRAPRepresentation(H{1},SumMatrixList(H(2:end)),prec);

end