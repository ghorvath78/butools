%  pi = DRPSolve(Q, prec)
%  
%  Computes the stationary solution of a discrete time 
%  Markov chain.
%  
%  Parameters
%  ----------
%  P : matrix, shape (M,M)
%      The matrix parameter of the rational process
%  prec : double, optional
%      Numerical precision for checking the rowsums.
%      The default value is 1e-14.
%      
%  Returns
%  -------
%  pi : row vector, shape (1,M)
%      The vector that satisfies 
%      `\pi\, P = \pi, \sum_i \pi_i=1`
%  
%  Notes
%  -----
%  Discrete time rational processes are like discrete time 
%  Markov chains, but the P matrix does not have to pass 
%  the :func:`CheckProbMatrix` test (but the rowsums still 
%  have to be ones).

function pi = DRPSolve (P)

    pi = CRPSolve(P-eye(size(P)));
end
