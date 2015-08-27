%  M = QBDFundamentalMatrices(B, L, F, matrices, precision, maxNumIt, method)
%  
%  Returns the fundamental matrices corresponding to the
%  given QBD. Matrices R, G and U can be returned
%  depending on the "matrices" parameter.
%  
%  The implementation is based on [1]_, please cite it if
%  you use this method.
%  
%  Parameters
%  ----------
%  B : matrix, shape (N,N)
%      The matrix corresponding to backward transitions
%  L : matrix, shape (N,N)
%      The matrix corresponding to local transitions
%  F : matrix, shape (N,N)
%      The matrix corresponding to forward transitions
%  matrices : string
%      Specifies which matrices are required. 'R' means 
%      that only matrix 'R' is needed. 'UG' means that
%      matrices U and G are needed. Any combinations of
%      'R', 'G' and 'U' are allowed, in any order.
%  precision : double, optional
%      The matrices are computed iteratively up to this
%      precision. The default value is 1e-14
%  maxNumIt : int, optional
%      The maximal number of iterations. The default value
%      is 50.
%  method : {"CR", "LR", "NI", "FI", "IS"}, optional
%      The method used to solve the matrix-quadratic
%      equation (CR: cyclic reduction, LR: logarithmic
%      reduction, NI: Newton iteration, FI: functional
%      iteration, IS: invariant subspace method). The 
%      default is "CR". "CR" is the only supported 
%      method in the Mathematica and Python implementation.
%  
%  Returns
%  -------
%  M : list of matrices
%      The list of calculated matrices in the order as
%      requested in the 'matrices' parameter.
%  
%  Notes
%  -----
%  Discrete and continuous QBDs are both supported, the
%  procedure auto-detects it based on the diagonal entries
%  of matrix L.
%  
%  References
%  ----------
%  .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
%         B. (2006, October). Structured Markov chains 
%         solver: software tools. In Proceeding from the
%         2006 workshop on Tools for solving structured 
%         Markov chains (p. 14). ACM.

function varargout = QBDFundamentalMatrices (B, L, F, matrices, precision, maxNumIt, method)

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
        fun = @QBD_CR;
    elseif strcmp(method,'LR')
        fun = @QBD_LR;
    elseif strcmp(method,'NI')
        fun = @QBD_NI;
    elseif strcmp(method,'IS')
        fun = @QBD_IS;
    elseif strcmp(method,'FI')
        fun = @QBD_FI;
    end        
    
    global BuToolsVerbose;
    
    if strfind(matrices,'R') 
        [G,R] = feval (fun, B, L, F, 'MaxNumIt', maxNumIt, 'Verbose', BuToolsVerbose);
    else
        G = feval (fun, B, L, F, 'MaxNumIt', maxNumIt, 'Verbose', BuToolsVerbose);
    end

    ret = cell(1,length(matrices));
    for i=1:length(matrices)
        if matrices(i)=='G'
            ret{i} = G;
        elseif matrices(i)=='R'
            ret{i} = R;
        elseif matrices(i)=='U'
            % we did not let smctools calculate U since there is a bug:
            % QBD_EG does not convert matrix U back to continuous 
            ret{i} = L+F*G;
        end
    end
    varargout=ret;
end
