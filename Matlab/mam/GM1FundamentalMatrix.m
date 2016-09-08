%  R = GM1FundamentalMatrix(A, precision, maxNumIt, method)
%  
%  Returns matrix R corresponding to the G/M/1 type Markov
%  chain given by matrices A.
%  
%  Matrix R is the minimal non-negative solution of the 
%  following matrix equation:
%  
%  .. math::
%      R = A_0 + R A_1 + R^2 A_2 + R^3 A_3 + \dots.
%  
%  The implementation is based on [1]_, please cite it if
%  you use this method.
%  
%  Parameters
%  ----------
%  A : length(M) list of matrices of shape (N,N)
%      Matrix blocks of the G/M/1 type generator in the 
%      regular part, from 0 to M-1.
%  precision : double, optional
%      Matrix R is computed iteratively up to this
%      precision. The default value is 1e-14
%  maxNumIt : int, optional
%      The maximal number of iterations. The default value
%      is 50.
%  method : {"CR", "RR", "NI", "FI", "IS"}, optional
%      The method used to solve the matrix-quadratic
%      equation (CR: cyclic reduction, RR: Ramaswami
%      reduction, NI: Newton iteration, FI: functional
%      iteration, IS: invariant subspace method). The 
%      default is "CR".
%  
%  Returns
%  -------
%  R : matrix, shape (N,N)
%      The R matrix of the G/M/1 type Markov chain.
%  
%  References
%  ----------
%  .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
%         B. (2006, October). Structured Markov chains 
%         solver: software tools. In Proceeding from the
%         2006 workshop on Tools for solving structured 
%         Markov chains (p. 14). ACM.

function R = GM1FundamentalMatrix (A, precision, maxNumIt, method)

    if ~exist('precision','var')
        precision = 1e-14;
    end
    
    if ~exist('maxNumIt','var')
        maxNumIt = 50;
    end
    
    if ~exist('method','var')
        method = 'CR';
    end   
    
    global BuToolsVerbose;

    N = size(A{1},2);
    Am = zeros(N, N*length(A));
    for i=1:length(A)
        Am(:,(i-1)*N+1:i*N) = A{i};
    end

    R = GIM1_R (Am, 'A', method, 'MaxNumIt', maxNumIt, 'Verbose', BuToolsVerbose, 'EpsilonValue', precision);
end
