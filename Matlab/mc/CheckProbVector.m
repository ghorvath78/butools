%  r = CheckProbVector(pi, sub, prec)
%  
%  Checks if the vector is a valid probability vector: the 
%  vector has only non-negative elements, the sum of the 
%  vector elements is 1.
%  
%  If parameter "sub" is set to true, it checks if the 
%  vector is a valid substochastic vector: the vector has 
%  only non-negative elements, the sum of the elements are
%  less than 1.
%  
%  Parameters
%  ----------
%  pi : vector, shape (1, M) or (M, 1)
%      The matrix to check.
%  sub : bool, optional
%      If false, the procedure checks for stochastic, if 
%      true, it checks for sub-stochastic property. The 
%      default value is false.
%  prec : double, optional
%      Numerical precision. Entries with absolute value 
%      less than prec are considered to be zeros. The 
%      default value is 1e-14.
%      
%  Returns
%  -------
%  r : bool
%      The result of the check.

function r = CheckProbVector (pi,sub,prec)

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

    if ~exist('sub','var')
        sub = 0;
    end

    r = false;
    
    if min(pi)<-prec
        if BuToolsVerbose
            fprintf ('CheckProbVector: The vector has negative element (precision: %g)!\n', prec);
        end
        return;
    end

    if sub
        if sum(pi)>1+prec*numel(pi)
            if BuToolsVerbose
                fprintf ('CheckProbVector: The sum of the substochastic vector is not less than 1 (precision: %g)!\n', prec);
            end
            return;
        end
    else
        if abs(sum(pi)-1)>prec*numel(pi)
            if BuToolsVerbose
                fprintf ('CheckProbVector: The sum of the vector is not 1 (precision: %g)!\n', prec);
            end
            return;
        end
    end
    
    r = true;
end