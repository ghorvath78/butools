%  pi = DTMCSolve(Q, prec)
%  
%  Computes the stationary solution of a discrete time 
%  Markov chain.
%  
%  Parameters
%  ----------
%  P : matrix, shape (M,M)
%      The transition probability matrix of the Markov 
%      chain
%  prec : double, optional
%      Numerical precision for checking whether P is a 
%      valid generator. The default value is 1e-14.
%      
%  Returns
%  -------
%  pi : row vector, shape (1,M)
%      The vector that satisfies `\pi\, P = \pi, \sum_i \pi_i=1`
%  
%  Notes
%  -----
%  The procedure raises an exception if :code:`butools.checkInput` 
%  is set to :code:`true` and :func:`CheckProbMatrix(P)` fails.

function pi = DTMCSolve (P)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end

    if BuToolsCheckInput && ~CheckProbMatrix(P, false)
        error('DTMCSolve: The given matrix is not a valid transient probability matrix. If you are sure you want this use DRPSolve instead of DTMCSolve.');
    end

    pi = DRPSolve(P);
end
