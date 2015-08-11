%  r = CheckMERepresentation(alpha, A, prec)
%  
%  Checks if the given vector and matrix define a valid matrix-
%  exponential representation.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      Initial vector of the matrix-exponential distribution 
%      to check
%  A : matrix, shape (M,M)
%      Matrix parameter of the matrix-exponential distribution
%      to check
%  prec : double, optional
%      Numerical precision. The default value is 1e-14.
%  
%  Returns
%  -------
%  r : bool
%      True, if the matrix is a square matrix, the vector and 
%      the matrix have the same size, the dominant eigenvalue
%      is negative and real
%  
%  Notes
%  -----
%  This procedure does not check the positivity of the density!
%  Call 'CheckMEPositiveDensity' if it is needed, but keep in
%  mind that it can be time-consuming, while this procedure
%  is fast.

function r = CheckMERepresentation (alpha, A, prec)

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if size(A,1)~=size(A,2)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The matrix is not a square matrix!\n');
        end
        r = false;
        return;
    end

    if length(alpha)~=size(A,1)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The vector and the matrix have different sizes!\n');
        end
        r = false;
        return;
    end

    if sum(alpha)<-prec*length(alpha) || sum(alpha)>1+prec*length(alpha)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The sum of the vector elements is less than zero or greater than one (precision: %g)!\n',prec);
        end
        r = false;
        return;
    end

    if max(real(eig(A)))>=prec
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: There is an eigenvalue of the matrix with non-negative real part (at precision %g)!\n',prec);
        end
        r = false;
        return;
    end
    
    ev = eig(A);
    [~,ix] = sort(abs(real(ev)));
    maxev = ev(ix(1));

    if ~isreal(maxev)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The dominant eigenvalue of the matrix is not real!\n');
        end
        r = false;
        return;
    end

    if sum(abs(ev(1:end))==abs(maxev)) > 1 && BuToolsVerbose
        fprintf ('CheckMERepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!\n');
    end

    r = true;
end
