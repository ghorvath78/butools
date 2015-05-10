%  r = CheckDMRAPRepresentation(H, prec)
%  
%  Checks if the input matrixes define a discrete time MRAP.
%  
%  All matrices H0...HK must have the same size, the dominant
%  eigenvalue of H0 is real and less than 1, and the rowsum of 
%  H0+H1+...+HK is 1 (up to the numerical precision).
%  
%  Parameters
%  ----------
%  H : list/cell of matrices, length(K)
%      The H0...HK matrices of the DMRAP to check
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckDMRAPRepresentation(H,prec)

    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    r = CheckDRAPRepresentation(H{1},SumMatrixList(H(2:end)),prec);

end