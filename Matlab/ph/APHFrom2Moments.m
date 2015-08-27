%  [alpha, A] = APHFrom2Moments(moms, maxSize)
%  
%  Returns an acyclic PH which has the same 2 moments as
%  given. If detects the order and the structure 
%  automatically to match the given moments.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(2)
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

function [alpha, A] = APHFrom2Moments (moms)
    cv2 = moms(2)/moms(1)^2 - 1.0;
    lambda = 1.0 / moms(1);
    N = max(ceil(1.0/cv2), 2);
    p = 1.0 / (cv2 + 1.0 + (cv2-1.0)/(N-1));
    A = -lambda*p*N * eye(N);
    for i=1:N-1
        A(i,i+1) = -A(i,i);
    end
    A(N,N) = -lambda*N;
    alpha = zeros(1,N);
    alpha(1) = p;
    alpha(N) = 1.0 - p;
end

