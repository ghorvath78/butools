%  r = CheckMGRepresentation(alpha, A, prec)
%  
%  Checks if the given vector and matrix define a valid matrix-
%  geometric representation.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      Initial vector of the matrix-geometric distribution 
%      to check
%  A : matrix, shape (M,M)
%      Matrix parameter of the matrix-geometric distribution
%      to check
%  prec : double, optional
%      Numerical precision. The default value is 1e-14.
%  
%  Returns
%  -------
%  r : bool
%      True, if the matrix is a square matrix, the vector and 
%      the matrix have the same size, the dominant eigenvalue
%      is positive, less than 1 and real. 
%  
%  Notes
%  -----
%  This procedure does not check the positivity of the density!
%  The discrete counterpart of 'CheckMEPositiveDensity' does
%  not exist yet (research is needed).

function r = CheckMGRepresentation(alpha, A, prec)

    global BuToolsVerbose
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if size(A,1) ~= size(A,2)
        if BuToolsVerbose
            fprintf('CheckMGRepresentation: The matrix is not a quadratic matrix!\n');
        end
        r = false;
        return;
    end

    if size(alpha,2) ~= size(A,1)
        if BuToolsVerbose
            fprintf('CheckMGRepresentation: The vector and the matrix have different sizes!\n');
        end
        r = false;
        return;
    end

    if sum(alpha)<-prec || sum(alpha)>1+prec
        if BuToolsVerbose
            fprintf ('CheckMGRepresentation: The sum of the vector elements is less than zero or greater than one (precision: %g)!\n',prec);
        end
        r = false;
        return;
    end

    ev = EigSort(eig(A));
    maxev = ev(1);

    if ~isreal(maxev)
        if BuToolsVerbose
            fprintf('CheckMGRepresentation: The largest eigenvalue of the matrix is complex!\n');
        end
        r = false;
        return;
    end

    if maxev > 1+prec
        if BuToolsVerbose
            fprintf('CheckMGRepresentation: The largest eigenvalue of the matrix is greater than 1 (precision: %g)!\n',prec);
        end
        r = false;
        return;
    end

    if sum(abs(ev(1:end))==abs(maxev)) > 1
         if BuToolsVerbose
            fprintf('CheckMGRepresentation Warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!\n');
         end
    end

    r = true;
end
