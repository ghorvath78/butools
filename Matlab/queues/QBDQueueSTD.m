%  [alpha, A] = QBDQueueSTD(B, L, F, L0, transToPH, prec)
%  
%  Returns the matrix-exponential distribution of the sojourn
%  time of the jobs in a QBD queue.
%  
%  QBD queues have a background continuous time Markov chain
%  with generator Q whose the transitions can be partitioned
%  into three sets: transitions accompanied by an arrival
%  of a new job (F, forward), transitions accompanied by 
%  the service of the current job in the server (B, 
%  backward) and internal transitions (L, local). 
%  Thus we have Q=B+L+F.
%  
%  Parameters
%  ----------
%  B : matrix, shape(N,N)
%      Transitions of the background process accompanied by 
%      the service of the current job in the server
%  L : matrix, shape(N,N)
%      Internal transitions of the background process 
%      that do not generate neither arrival nor service
%  F : matrix, shape(N,N)
%      Transitions of the background process accompanied by 
%      an arrival of a new job
%  L0 : matrix, shape(N,N)
%      Internal transitions of the background process when
%      there are no jobs in the queue
%  transToPH : bool, optional
%      If true, the result is transformed to a phase-type
%      representation. The default value is false
%  prec : double, optional
%      Numerical precision to check if the input is valid and
%      it is also used as a stopping condition when solving
%      the matrix-quadratic equation
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-exponential 
%      distribution representing the sojourn time of the jobs
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      distribution representing the sojourn time of the jobs
%  
%  Notes
%  -----
%  While the transformation of the results to phase-type
%  representation is always possible theoretically, it may
%  introduce numerical problems in some cases.

function [alpha, A] = QBDQueueSTD (B, L, F, L0, transToPH, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('transToPH','var')
        transToPH = false;
    end

    % obtain matrix-geometric representation
    [pi0, R] = QBDSolve (B, L, F, L0, prec);

    N = length(pi0);
    I = eye(N);
    
    U = L + R*B;
    Rh = inv(-U)*F;
    eta = pi0*F*inv(I-Rh);
    eta = eta/sum(eta);
    z = reshape(I,N*N,1);
   if transToPH
        % transform to ph distribution
        ix = (1:N);
        nz = ix(eta>prec);
        Delta = diag(eta);
        A = kron(L+F,I(nz,nz)) + kron(B,inv(Delta(nz,nz))*Rh(nz,nz)'*Delta(nz,nz));
        alpha = z'*kron(I,Delta(:,nz));
   else
        % transform it such that the closing vector is a vector of ones
        % this is the way butools accepts ME distributions
        Bm = TransformToOnes(z);
        Bmi = inv(Bm);
        A = Bm * (kron(L'+F',I) + kron(B',Rh)) * Bmi;
        alpha = kron(ones(1,N), eta) * Bmi;
   end
end
