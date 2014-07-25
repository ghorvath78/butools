%  [alpha, A] = MAPMAP1QLD(D0, D1, S0, S1, transToPH, prec)
%  
%  Returns the matrix-geometric distribution of the queue
%  length of a MAP/MAP/1 queue.
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

function [alpha, A] = MAPMAP1QLD (D0, D1, S0, S1, transToPH, prec)

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
        error('MAPMAP1QLD: The arrival process (D0,D1) is not a valid MAP representation!');
    end
    
    if BuToolsCheckInput && ~CheckMAPRepresentation(S0,S1,prec)
        error('MAPMAP1QLD: The service process (S0,S1) is not a valid MAP representation!');
    end

    IA = eye(size(D0,1));
    IS = eye(size(S0,1));
    
    [alpha, A] = QBDQueueQLD (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), transToPH, prec);    
end

