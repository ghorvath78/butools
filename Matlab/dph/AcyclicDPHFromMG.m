%  [beta, B] = AcyclicDPHFromMG(alpha, A, precision)
%  
%  Transforms a matrix-geometric representation to an acyclic
%  DPH representation of the same size, if possible.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,N)
%      Initial vector of the distribution
%  A : matrix, shape (N,N)
%      Matrix parameter of the distribution
%  precision : double, optional
%      Vector and matrix entries smaller than the precision
%      are considered to be zeros. The default value is 1e-14.
%  
%  Returns
%  -------
%  beta : matrix, shape (1,M)
%      The initial probability vector of the acyclic discrete
%      phase-type representation
%  B : matrix, shape (M,M)
%      Transition probability matrix of the acyclic discrete
%      phase-type representation
%  
%  Notes
%  -----
%  Contrary to 'AcyclicPHFromME' of the 'ph' package, this 
%  procedure is not able to extend the size in order to obtain
%  a Markovian initial vector.
%  
%  Raises an error if A has complex eigenvalues. In this case
%  the transformation to an acyclic representation is not 
%  possible

function [beta,B] = AcyclicDPHFromMG (alpha, A, prec)

    if ~exist('prec','var')
        prec = 10^-14;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation(alpha, A)
        error('AcyclicDPHFromMG: Input isn''t a valid MG distribution!');
    end

    lambda = eig(A);

    if max(abs(imag(lambda)))>prec
        error('AcyclicDPHFromMG: The input matrix has complex eigenvalue!');
    end

    lambda2 = EigSort(lambda);
    lambda3 = lambda2(lambda2~=lambda2(end));
    mx = diag(lambda2)+diag(1-lambda3, 1);
    T = SimilarityMatrix (A, mx);
    beta = alpha*T;
    B = mx;

    if ~CheckDPHRepresentation (beta, B, prec)
        error('AcyclicDPHFromMG: No acyclic representation found!');
    end
end
