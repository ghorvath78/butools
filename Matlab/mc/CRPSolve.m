%  pi = CRPSolve(Q, prec)
%  
%  Computes the stationary solution of a continuous time 
%  rational process (CRP).
%  
%  Parameters
%  ----------
%  Q : matrix, shape (M,M)
%      The generator matrix of the rational process
%  prec : double, optional
%      Numerical precision for checking the rowsums.
%      The default value is 1e-14.
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

function pi = CRPSolve (Q, prec)

    if ~exist('prec','var')
        prec=1e-14;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end

    if BuToolsCheckInput && any(sum(Q,2)>prec)
        error('CRPSolve: The matrix has a rowsum which isn''t zero!');
    end

    M = Q;
    M(:,1) = 1;
    m = zeros(1,size(M,1));
    m(1) = 1.0;
    pi = m/M;
end
