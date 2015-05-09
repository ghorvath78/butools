%  [beta, B] = AcyclicPHFromME(alpha, A, maxSize, precision)
%  
%  Transforms an arbitrary matrix-exponential representation
%  to an acyclic phase-type representation. (see [1]_).
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,N)
%      Initial vector of the distribution
%  A : matrix, shape (N,N)
%      Matrix parameter of the distribution
%  maxSize : int, optional
%      The maximum number of phases for the result.
%      The default value is 100.
%  precision : double, optional
%      Vector and matrix entries smaller than the precision
%      are considered to be zeros. The default value is 1e-14.
%  
%  Returns
%  -------
%  beta : matrix, shape (1,M)
%      The initial probability vector of the Markovian 
%      acyclic representation
%  B : matrix, shape (M,M)
%      Transient generator matrix of the Markovian 
%      acyclic representation
%  
%  Notes
%  -----
%  Raises an error if no Markovian acyclic representation
%  has been found.
%  
%  References
%  ----------
%  .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
%          phase-type distributions," Stoch. Models 15, 759-778 
%          (1999)

function [beta, B] = AcyclicPHFromME (alpha, A, maxSize, precision)
        
    if ~exist('precision','var')
        precision = 1e-14;
    end

    if ~exist('maxSize','var')
        maxSize = 100;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('AcyclicPHFromME: Input isn''t a valid ME distribution!');
    end

    G = TransformToAcyclic (A, maxSize, precision);

    % find transformation matrix
    T = SimilarityMatrix (A, G);
    gamma = real(alpha*T);

    if min(gamma) > -precision
        beta = gamma;
        B = G;
    else
        [beta, B] = ExtendToMarkovian (gamma, G, maxSize, precision);
    end
    
    if ~CheckPHRepresentation(beta, B, precision)
        error('AcyclicPHFromME: No acyclic representation found up to the given size and precision!');
    end   
end
