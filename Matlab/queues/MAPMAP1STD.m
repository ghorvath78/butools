%  [alpha, A] = MAPMAP1STD(D0, D1, S0, S1, transToPH, prec)
%  
%  Returns the matrix-exponential distribution of the sojourn
%  time of the jobs in a MAP/MAP/1 queue.
%  
%  In a MAP/MAP/1 queue both the arrival and the service
%  processes are characterized by Markovian arrival 
%  processes.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape(N,N)
%      The transitions of the arrival MAP not accompanied by
%      job arrivals
%  D1 : matrix, shape(N,N)
%      The transitions of the arrival MAP accompanied by
%      job arrivals
%  S0 : matrix, shape(N,N)
%      The transitions of the service MAP not accompanied by
%      job service
%  S1 : matrix, shape(N,N)
%      The transitions of the service MAP accompanied by
%      job service
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
%  
%  Notes
%  -----
%  While the transformation of the results to phase-type
%  representation is always possible theoretically, it may
%  introduce numerical problems in some cases.

function [alpha, A] = MAPMAP1STD (D0, D1, S0, S1, transToPH, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end
    
    if ~exist('transToPH','var')
        transToPH = false;
    end
    
    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1,prec)
        error('MAPMAP1STD: The arrival process (D0,D1) is not a valid MAP representation!');
    end
    
    if BuToolsCheckInput && ~CheckMAPRepresentation(S0,S1,prec)
        error('MAPMAP1STD: The service process (S0,S1) is not a valid MAP representation!');
    end

    IA = eye(size(D0,1));
    IS = eye(size(S0,1));
   
    B = kron(IA,S1);
    L = kron(D0,IS)+kron(IA,S0);
    F = kron(D1,IS);
    
    [pi0, R] = QBDSolve (B, L, F, kron(D0,IS), prec);
    
    N = length(pi0);
    I = eye(N);

    % calculate the distribution of the age at departures
    U = L + R*B;
    Rh = inv(-U)*F;
    T = kron(IA,S0) + Rh * B;
    eta = pi0*F*inv(I-Rh);
    eta = eta/sum(eta);

    if transToPH
        % transform it to PH representation
        beta = CTMCSolve(S0+S1);
        theta = DTMCSolve(inv(-D0)*D1);
        vv = kron(theta,beta);
        ix = 1:N;
        nz = ix(vv>prec);
        delta = diag(vv(nz));   
        alpha = ones(1,N)*B(nz,:)'*delta / sum(beta*S1);
        A = inv(delta)*T(nz,nz)'*delta;
    else
        % transform it to ME
        A = T;
        alpha = eta;
    end
end

