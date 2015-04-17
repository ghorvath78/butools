%  [pi0, R] = QBDSolve (B, L, F, L0, prec)
%  
%  Returns the parameters of the matrix-geometrically 
%  distributed stationary distribution of a QBD.
%  
%  Using vector pi0 and matrix R provided by this function
%  the stationary solution can be obtained by
%  
%  .. math::
%      \pi_k=\pi_0 R^k.    
%  
%  Parameters
%  ----------
%  B : matrix, shape (N,N)
%      The matrix corresponding to backward transitions
%  L : matrix, shape (N,N)
%      The matrix corresponding to local transitions
%  F : matrix, shape (N,N)
%      The matrix corresponding to forward transitions
%  L0 : matrix, shape (N,N)
%      The matrix corresponding to local transitions at
%      level zero
%  precision : double, optional
%      The fundamental matrix R is computed up to this
%      precision. The default value is 1e-14
%  
%  Returns
%  -------
%  pi0 : matrix, shape (1,N)
%      The stationary probability vector of level zero
%  R : matrix, shape (N,N)
%      The matrix parameter of the matrix geometrical
%      distribution of the QBD 

function [pi0, R] = QBDSolve (B, L, F, L0, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    m = size(L0,1);
    I = eye(m);

    R = QBDFundamentalMatrices (B, L, F, 'R', prec);
    
    % Convert to discrete time problem, if needed
    if sum(diag(L0)) < 0 % continues time
        lamb = max(-diag(L0));
        B = B / lamb;
        L0 = L0 / lamb + I;
    end
    
    pi0 = DTMCSolve(L0+R*B);
    nr = sum(pi0*inv(I-R));
    pi0 = pi0 / nr;
end
