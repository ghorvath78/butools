%  pi = QBDStationaryDistr (pi0, R, K)
%  
%  Returns the stationary distribution of a QBD up to a
%  given level K.
%  
%  Parameters
%  ----------
%  pi0 : matrix, shape (1,N)
%      The stationary probability vector of level zero
%  R : matrix, shape (N,N)
%      The matrix parameter of the matrix geometrical
%      distribution of the QBD 
%  K : integer
%      The stationary distribution is returned up to
%      this level.
%  
%  Returns
%  -------
%  pi : array, length (K+1)*N
%      The stationary probability vector up to level K

function pi = QBDStationaryDistr (pi0, R, K)

    m = size(R,1);
    pi = zeros(1,(K+1)*m);
    pi(1:m) = pi0;
    pix = pi0;
    for k=1:K
        pix = pix * R;
        pi(k*m+1:(k+1)*m) = pix;
    end
end
