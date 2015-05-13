%  pi = CRPSolve(Q)
%  
%  Computes the stationary solution of a continuous time 
%  rational process (CRP).
%  
%  Parameters
%  ----------
%  Q : matrix, shape (M,M)
%      The generator matrix of the rational process
%      
%  Returns
%  -------
%  pi : row vector, shape (1,M)
%      The vector that satisfies 
%      `\pi\, Q = 0, \sum_i \pi_i=1`
%  
%  Notes
%  -----
%  Continuous time rational processes are like continuous 
%  time Markov chains, but the generator does not have to 
%  pass the :func:`CheckGenerator` test (but the rowsums 
%  still have to be zeros).

function pi = CRPSolve (Q)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end

    if BuToolsCheckInput && any(abs(sum(Q,2))>BuToolsCheckPrecision)
        error('CRPSolve: The matrix has a rowsum which isn''t zero!');
    end

    M = Q;
    M(:,1) = 1;
    m = zeros(1,size(M,1));
    m(1) = 1.0;
    pi = m/M;
end
