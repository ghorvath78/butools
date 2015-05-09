%  [alpha, A] = APHFrom3Moments(moms, maxSize)
%  
%  Returns an acyclic PH which has the same 3 moments as
%  given. If detects the order and the structure 
%  automatically to match the given moments.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(3)
%    The moments to match
%  maxSize : int, optional
%    The maximal size of the resulting APH. The default value
%    is 100.
%  
%  Returns
%  -------
%  alpha : vector, shape (1,M)
%    The initial probability vector of the APH
%  A : matrix, shape (M,M)
%    Transient generator matrix of the APH
%  
%  Raises an error if the moments are not feasible with an
%  APH of size "maxSize".
%  
%  References
%  ----------
%  .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
%         moments with minimal acyclic phase type 
%         distributions," Stochastic models, pp. 303-326, 2005.

function [alpha, A] = APHFrom3Moments (moms, maxSize, prec)

    if ~exist('prec','var')
        prec=1e-14;
    end
    if ~exist('maxSize','var')
        maxSize = 100;
    end
    
    
    m1 = moms(1);
    m2 = moms(2);
    m3 = moms(3);
    
    % detect number of phases needed
    n = 2;
    while n<maxSize && (APH2ndMomentLowerBound(m1, n) > m2 || APH3rdMomentLowerBound(m1, m2, n) >= m3 || APH3rdMomentUpperBound(m1, m2, n) <= m3)
        n = n + 1;
    end

    % if PH is too large, adjust moment to bounds
    if APH2ndMomentLowerBound(m1, n) > m2
        m2 = APH2ndMomentLowerBound(m1, n);
    end

    if APH3rdMomentLowerBound(m1, m2, n) > m3
        m3 = APH3rdMomentLowerBound(m1, m2, n);
    end

    if APH3rdMomentUpperBound(m1, m2, n) < m3
        m3 = APH3rdMomentUpperBound(m1, m2, n);
    end

    % compute normalized moments
    nmoms = NormMomsFromMoms ([m1, m2, m3]);
    n1 = nmoms(1);
    n2 = nmoms(2);
    n3 = nmoms(3);


    if n2>2 || n3 < 2*n2 - 1
        b = real(2*(4-n*(3*n2-4)) / (n2*(4+n-n*n3) + sqrt(n*n2)*sqrt(12*n2^2*(n+1)+16*n3*(n+1)+n2*(n*(n3-15)*(n3+1)-8*(n3+3)))));
        a = (b*n2-2)*(n-1)*b / (b-1) / n;
        p = (b-1) / a;
        lambda = (p*a+1) / n1;
        mu = (n-1)*lambda / a;
        % construct representation
        alpha = zeros(1,n);
        alpha(1) = p;
        alpha(n) = 1.0-p;
        A = zeros(n,n);
        A(n,n) = -lambda;
        for i=1:n-1
            A(i,i) = -mu;
            A(i,i+1) = mu;
        end
        return;
    else
        c4 = n2*(3*n2-2*n3)*(n-1)^2;
        c3 = 2*n2*(n3-3)*(n-1)^2;
        c2 = 6*(n-1)*(n-n2);
        c1 = 4*n*(2-n);
        c0 = n*(n-2);
        fs = roots([c4 c3 c2 c1 c0]);
        for i=1:length(fs)
            f = fs(i);
            if abs((n-1)*(n2*f^2-2*f+2)-n)<prec
                continue;
            end
            a = 2*(f-1)*(n-1) / ((n-1)*(n2*f^2-2*f+2)-n);
            p = (f-1)*a;
            lambda = (a+p) / n1;
            mu = (n-1) / (n1 - p/lambda);
            if isreal(p) && isreal(lambda) && isreal(mu)&& p>=0 && p<=1 && lambda>0 && mu>0
                alpha = zeros(1,n);
                alpha(1) = p;
                alpha(2) = 1-p;
                A = zeros(n,n);
                A(1,1) = -lambda;
                A(1,2) = lambda;
                for j=2:n
                    A(j,j) = -mu;
                    if j<n
                        A(j,j+1) = mu;
                    end
                end
                return;
            end
        end
    end
    error('No APH found for the given 3 moments!');
end
