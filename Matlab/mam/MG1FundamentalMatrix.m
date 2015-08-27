%  G = MG1FundamentalMatrix(A, precision, maxNumIt, method)
%  
%  Returns matrix G corresponding to the M/G/1 type Markov
%  chain defined by matrices A.
%  
%  Matrix G is the minimal non-negative solution of the 
%  following matrix equation:
%  
%  .. math::
%      G = A_0 + A_1 G + A_2 G^2 + A_3 G^3 + \dots.
%  
%  The implementation is based on [1]_, please cite it if
%  you use this method.
%  
%  Parameters
%  ----------
%  A : length(M) list of matrices of shape (N,N)
%      Matrix blocks of the M/G/1 type generator from
%      0 to M-1.
%  precision : double, optional
%      Matrix G is computed iteratively up to this
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
%  G : matrix, shape (N,N)
%      The G matrix of the M/G/1 type Markov chain.
%      (G is stochastic.)
%  
%  References
%  ----------
%  .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
%         B. (2006, October). Structured Markov chains 
%         solver: software tools. In Proceeding from the
%         2006 workshop on Tools for solving structured 
%         Markov chains (p. 14). ACM.

function G = MG1FundamentalMatrix (A, precision, maxNumIt, method)

    if ~exist('precision','var')
        precision = 1e-14;
    end
    
    if ~exist('maxNumIt','var')
        maxNumIt = 50;
    end
    
    if ~exist('method','var')
        method = 'CR';
    end
    
    if strcmp(method,'CR')
        fun = @MG1_CR;
    elseif strcmp(method,'RR')
        fun = @MG1_RR;
    elseif strcmp(method,'NI')
        fun = @MG1_NI;
    elseif strcmp(method,'IS')
        fun = @MG1_IS;
    elseif strcmp(method,'FI')
        fun = @MG1_FI;
    end        
    
    global BuToolsVerbose;
    
    N = size(A{1},2);
    Am = zeros(N, N*length(A));
    for i=1:length(A)
        Am(:,(i-1)*N+1:i*N) = A{i};
    end
    
    G = feval (fun, Am, 'MaxNumIt', maxNumIt, 'Verbose', double(BuToolsVerbose), 'EpsilonValue', precision);
end
