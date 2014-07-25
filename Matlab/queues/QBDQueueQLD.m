%  [alpha, A] = QBDQueueQLD(B, L, F, L0, transToPH, prec)
%  
%  Returns the matrix-geometric distribution of the queue
%  length of a QBD queue.
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
%      If true, the result is transformed to a discrete 
%      phase-type representation. The default value is false
%  prec : double, optional
%      Numerical precision to check if the input is valid and
%      it is also used as a stopping condition when solving
%      the matrix-quadratic equation
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-geometric 
%      distribution of the number of jobs in the system
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-geometric
%      distribution of the number of jobs in the system
%  
%  Notes
%  -----
%  While the transformation of the results to a discrete 
%  phase-type representation is always possible 
%  theoretically, it may introduce numerical problems in 
%  some cases.

function [alpha, A] = QBDQueueQLD (B, L, F, L0, transToPH, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('transToPH','var')
        transToPH = false;
    end

    % obtain matrix-geometric representation
    [pi0, R] = QBDSolve (B, L, F, L0, prec);

    N = length(pi0);
    if transToPH
        % transform it to DPH
        alpha = pi0*R*inv(eye(N)-R);
        A = inv(diag(alpha))*R'*diag(alpha);
    else
        % transform it to MG
        B = TransformToOnes(sum(inv(eye(N)-R)*R,2));
        Bi = inv(B);
        A = B*R*Bi;
        alpha = pi0*Bi;       
    end
end

