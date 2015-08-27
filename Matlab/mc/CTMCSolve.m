%  pi = CTMCSolve(Q)
%  
%  Computes the stationary solution of a continuous time 
%  Markov chain.
%  
%  Parameters
%  ----------
%  Q : matrix, shape (M,M)
%      The generator matrix of the Markov chain
%      
%  Returns
%  -------
%  pi : row vector, shape (1,M)
%      The vector that satisfies `\pi\, Q = 0, \sum_i \pi_i=1`
%  
%  Notes
%  -----
%  The procedure raises an exception if :code:`checkInput` 
%  is set to :code:`true` and :func:`CheckGenerator` (Q) fails.

function pi = CTMCSolve(Q)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end

    if BuToolsCheckInput && ~CheckGenerator(Q, false)
        error('CTMCSolve: The given matrix is not a valid generator. If you are sure you want this use CRPSolve instead of CTMCSolve.');
    end

    pi = CRPSolve (Q);
end

