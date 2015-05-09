%  [beta, B] = MinimalRepFromME(alpha, A, how, precision)
%  
%  Returns the minimal representation of the given ME 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-exponential 
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      distribution.
%  how : {"obs", "cont", "obscont", "moment"}, optional        
%      Determines how the representation is minimized. 
%      Possibilities:
%      'obs': observability, 
%      'cont': controllability,
%      'obscont': the minimum of observability and 
%      controllability order,
%      'moment': moment order (which is the default).
%  precision : double, optional
%     Precision used by the Staircase algorithm. The default
%     value is 1e-12.
%  
%  Returns
%  -------
%  beta : vector, shape (1,N)
%      The initial vector of the minimal representation
%  B : matrix, shape (N,N)
%      The matrix parameter of the minimal representation
%  
%  References
%  ----------
%  .. [1]  P. Buchholz, M. Telek, "On minimal representation
%          of rational arrival processes." Madrid Conference on
%          Qeueuing theory (MCQT), June 2010.

function [beta, B] = MinimalRepFromME (alpha, A, how, precision)

    if ~exist('precision','var')
        precision=1e-12;
    end
    if ~exist('how','var')
        how='moment';
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('MinimalRepFromME: Input isn''t a valid ME distribution!');
    end

    if strcmp(how,'cont')
        H0 = A;
        H1 = sum(-A,2) * alpha;
        [B, n] = MStaircase ({H0, H1}, ones(size(A,1),1), precision);
        beta = alpha*B;
        beta = beta(1:n);
        B = inv(B)*A*B;
        B = B(1:n,1:n);
    elseif strcmp(how,'obs')
        H0 = A;
        H1 = sum(-A,2) * alpha;        
        [B, n] = MStaircase ({H0',H1'}, alpha', precision);
        beta = alpha*B;
        beta = beta(1:n);
        B = inv(B)*A*B;
        B = B(1:n,1:n);
    elseif strcmp(how,'obscont')
        [alphav, Av] = MinimalRepFromME (alpha, A, 'cont', precision);
        [beta, B] = MinimalRepFromME (alphav, Av, 'obs', precision);
    elseif strcmp(how,'moment')
        mo = MEOrder (alpha, A, 'moment', precision);
        [beta, B] = MEFromMoments(MomentsFromME(alpha, A, 2*mo-1));
    else
        error('MinimalRepFromME: Invalid ''how'' parameter!')
    end
end

