function [alpha, A] = APHFrom2Moments (m1, m2)
    cv2 = m2/m1^2 - 1.0;
    lambda = 1.0 / m1;
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

