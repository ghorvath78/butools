%  r = CheckRAPRepresentation(H0, H1, prec)
%  
%  Checks if the input matrixes define a continuous time RAP.
%  
%  Matrices H0 and H1 must have the same size, the dominant
%  eigenvalue of H0 is negative and real, and the rowsum of 
%  H0+H1 is 0 (up to the numerical precision).
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the RAP to check
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the RAP to check
%  prec : double, optional
%      Numerical precision, the default value is 1e-14
%  
%  Returns
%  -------
%  r : bool 
%      The result of the check

function r = CheckRAPRepresentation (D0, D1, prec)

    global BuToolsVerbose;

    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if size(D0,1)~=size(D0,2)
         if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: D0 is not a quadratic matrix!\n');
        end
        r = false;
        return;
    end

    if size(D1,1)~=size(D1,2)
        if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: D1 is not a quadratic matrix!\n');
        end
        r = false;
        return;
    end

    if size(D0)~=size(D1)
        if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: D0 and D1 have different sizes!\n');
        end
        r = false;
        return;
    end

    if max(abs((D0+D1)*ones(size(D0,1),1))) > prec
        if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: A rowsum of D0+D1 is not 0!(precision: %g)\n',prec);
        end
        r = false;
        return;
    end


    ev=eig(D0);
    if max(real(ev))>=-prec
        if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part (at precision %g)\n',prec);
        end
        r = false;
        return;
    end

    ceig=ev;
    reig=ev;
    for i=1:length(ev)
        if isreal(ev(i))
            ceig(i)=-Inf;
            reig(i)=real(ev(i));
        else
            ceig(i)=real(ev(i));
            reig(i)=-Inf;
        end
    end

    if max(reig) < max(ceig)
        if BuToolsVerbose
            fprintf ('CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!\n');
        end
        r = false;
        return;
    end

    if max(reig)==max(ceig)
        if BuToolsVerbose
            fprintf('CheckRAPRepresentation: The dominant and a complex eigenvalue of D0 has the same real part!\n');
        end
    end

    r = true;
end
