%  [alpha, A] = MAPPH1STD(D0, D1, sigma, S, transToPH, prec)
%  
%  Returns the matrix-exponential distribution of the sojourn
%  time of the jobs in a MAP/PH/1 queue.
%  
%  In a MAP/PH/1 queue the arrival processes is given by a 
%  Markovian arrival processes, and the service time follows
%  a phase-type distribution.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape(N,N)
%      The transitions of the arrival MAP not accompanied by
%      job arrivals
%  D1 : matrix, shape(N,N)
%      The transitions of the arrival MAP accompanied by
%      job arrivals
%  sigma : matrix, shape(1,N)
%      The initial probability vector of the phase-type
%      distribution corresponging to the service time
%  S : matrix, shape(N,N)
%      The transient generator of the phase-type 
%      distribution corresponging to the service time
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

function [alpha, A] = MAPPH1STD (D0, D1, sigma, S, transToPH, prec)

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
        error('MAPPH1STD: The arrival process (D0,D1) is not a valid MAP representation!');
    end
    
    if BuToolsCheckInput && ~CheckPHRepresentation(sigma,S,prec)
        error('MAPPH1STD: The service distribution (sigma,S) is not a valid PH representation!');
    end

    [alpha, A] = MAPMAP1STD (D0, D1, S, sum(-S,2)*sigma, transToPH, prec);    
end

