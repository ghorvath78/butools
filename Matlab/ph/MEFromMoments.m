%  [alpha, A] = MEFromMoments(moms)
%  
%  Creates a matrix-exponential distribution that has the
%  same moments as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(2*M-1)
%      The list of moments. The order of the resulting 
%      matrix-exponential distribution is 
%      determined based on the number of moments given. To 
%      obtain a matrix exponential distribution of order M,
%      2*M-1 moments are required.
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-exponential 
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      distribution.
%  
%  References
%  ----------
%  .. [1] A. van de Liefvoort. The moment problem for 
%         continuous distributions. Technical report, 
%         University of Missouri, WP-CM-1990-02, Kansas City,
%         1990.

function [alpha, A] = MEFromMoments (moms)

    function K = appie (rmom)
        m = length (rmom);
        if rem(m,2)==0
            rm = rmom(1:m-1);
            m = m / 2;
        else
            rm = rmom;
            m = ceil(m/2);
        end
        rm = [1 rm];
        f = zeros(2*m, 1);
        f(1) = 1;
        y = zeros(2*m, 1);
        dd = zeros(2*m, 2*m);
        n = 0;
        k = 0;
        q = 1;
        d = zeros(m, 1);
        alpha = zeros(m, m);
        beta = zeros(m, 1);
        for i=2:2*m
            dd(i,i-1) = 1;
        end
        for i=1:2*m
            ro = q*rm*f;
            nold = n;
            n = nold + 1;
            yold = y;
            if n>0 && ro~=0
                if k>0
                    beta(k) = ro / (rm(2).^(d(k)+n-1));
                end
                k = k+1;
                d(k) = n;
                n = -n;
                q = q / ro;
                y = dd*f;
            elseif n<=0
                j = nold + d(k) + 1;
                alpha(k,j) = ro / rm(2).^(j-1);
            end
            f = dd*f - ro*yold;
        end

        if sum(d)~=m
            error ('Insufficient matrix order!');
        end

        K = zeros(m,m);
        K(1,1) = rm(2);
        for i=1:m-1
            K(i,i+1) = rm(2);
        end
        ind = d(1);
        for i=2:m
            if ind<m
                inc = d(i);
                ind = ind + inc;
                if ind<=m
                    K(ind, ind-inc-d(i-1)+1) = beta(i-1);
                    for j=1:inc
                        K(ind, ind-j+1) = alpha(i,j);
                    end
                end
            end
        end
    end

    K = appie (ReducedMomsFromMoms(moms));
    N = ceil(length(moms)/2);

    T = zeros(N,N);
    for i=1:N
        for j=1:i
            T(i,j) = 1;
        end
    end

    U = zeros(N,N);
    for i=1:N
        for j=i:N
            U(i,j) = 1 / (N-i+1);
        end
    end

    alpha = zeros(1,N);
    alpha(1) = 1;
    alpha = alpha * inv(T) * U;
    A = inv(-inv(U)*T*K*inv(T)*U);
end
