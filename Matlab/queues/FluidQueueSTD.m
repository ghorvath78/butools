%  [alpha, A] = FluidQueueSTD (Q, Rin, Rout, Q0, transToPH, prec)
%  
%  Returns the parameters of the matrix-exponentially 
%  distributed stationary sojourn time distribution of
%  fluid drops in a fluid queue.
%  
%  In a fluid queue there is a background continuous time
%  Markov chain (given by generator Q), and diagonal
%  matrix Rin (Rout) whose ith entry provides the 
%  fluid rate at which fluid enters the queue (can be 
%  served) while the background process is in state i.
%  
%  
%  Parameters
%  ----------
%  Q : matrix, shape (N,N)
%      The generator of the background Markov chain
%  Rin : matrix, shape (N,N)
%      Diagonal matrix containing the fluid input rates
%      associated to the states of the background process
%  Rout : matrix, shape (N,N)
%      Diagonal matrix containing the fluid output rates
%      associated to the states of the background process
%  Q0 : matrix, shape (N,N), optional
%      The generator of the background Markov chain when
%      the queue level is zero. If empty or missing, Q0=Q
%      is assumed
%  transToPH : bool, optional
%      If true, the result is transformed to a phase-type
%      representation, otherwise it is returned as a matrix
%      exponential representation. The default value is 
%      true
%  prec : double, optional
%      Numerical precision used to check wether the input 
%      is valid and it also serves as a stopping condition
%      when the algebraic Riccati equation is solved. The
%      delault value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-exponential 
%      sojourn time distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      sojourn time distribution.
%  
%  Notes
%  -----
%  While it is always possible to transform the result to 
%  a phase-type representation theoretically, this 
%  transformation step can be sensitive numerically.
%  
%  Notes
%  -----
%  The result can be very large! If the input and output
%  processes are independent, use the 'FluFluSTD' function
%  which much faster and returns a small representation.

function [alpha, A] = FluidQueueSTD (Q, Rin, Rout, Q0, transToPH, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('Q0','var')
        Q0 = [];
    end

    if ~exist('transToPH','var')
        transToPH = false;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckGenerator(Q,false,prec)
        error('FluidQueueSTD: Generator matrix Q is not Markovian!');
    end
    if BuToolsCheckInput && ~isempty(Q0) && ~CheckGenerator(Q0,false,prec)
        error('FluidQueueSTD: Generator matrix Q0 is not Markovian!');
    end
    if BuToolsCheckInput && (any(diag(Rin)<-prec) || any(diag(Rout)<-prec))
        error('FluidQueueSTD: Fluid rates Rin and Rout must be non-negative !');
    end  
    
    % obtain solution
    [mass0, ini, K, clo] = GeneralFluidSolve (Q, Rin-Rout, Q0, prec);
   
    N = size(Q,1);
    iniKi = linsolve(K',-ini')'; % iniki = ini*inv(-K);
    lambda = sum(mass0*Rin + iniKi*clo*Rin);
    if transToPH
        % transform it to PH
        Delta = diag(iniKi/lambda);
        alpha = reshape(clo*Rin,1,N*length(ini))*kron(eye(N),Delta);
        A = kron(Rout, inv(Delta)*K'*Delta) + kron(Q, eye(size(K)));
    else
        B = TransformToOnes(reshape(inv(-K)*clo*Rin,N*length(ini),1));
        Bi = inv(B);
        alpha = kron(ones(1,N), ini/lambda)*Bi;
        A = B*(kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K)))*Bi;        
    end
end

