%  order = MEOrder(alpha, A, kind, prec)
%  
%  Returns the order of the ME distribution (which is not 
%  necessarily equal to the size of the representation).
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-exponential 
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      distribution.
%  kind : {'obs', 'cont', 'obscont', 'moment'}, optional
%      Determines which order is computed. Possibilities: 
%      'obs': observability, 
%      'cont': controllability,
%      'obscont': the minimum of observability and 
%      controllability order,
%      'moment': moment order (which is the default).
%  prec : double, optional
%      Precision used to detect if the determinant of the 
%      Hankel matrix is zero (in case of kind="moment" only),
%      or the tolerance for the rank calculation. The
%      default value is 1e-10.
%  
%  Returns
%  -------
%  order : int
%      The order of ME distribution
%  
%  References
%  ----------
%  .. [1]  P. Buchholz, M. Telek, "On minimal representation
%          of rational arrival processes." Madrid Conference on
%          Qeueuing theory (MCQT), June 2010.

function order = MEOrder (alpha, A, kind, prec)

   if ~exist('prec','var')
       prec = 1e-10;
   end
   
   if ~exist('kind','var')
       kind = 'moment';
   end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('MinimalRepFromME: Input isn''t a valid ME distribution!');
    end
   
   N = length(alpha);  
   if strcmp(kind,'cont')
       re = zeros (N);
       for n=1:N
           re(n,:) = sum(A'^(n-1), 1);
       end
       order = rank (re, prec);
   elseif strcmp(kind,'obs')
       re = zeros (N);
       for n=1:N
           re(n,:) = alpha*A^(n-1);
       end
       order = rank (re, prec);
   elseif strcmp(kind,'obscont')
       re = zeros (N);
       for n=1:N
           re(n,:) = alpha*A^(n-1);
       end
       obsOrder = rank (re, prec);
       for n=1:N
           re(n,:) = sum(A'^(n-1), 1);
       end
       contOrder = rank (re, prec);
       order = min(obsOrder,contOrder);
   elseif strcmp(kind,'moment')
       order = MEOrderFromMoments (MomentsFromME (alpha, A), prec);
   else
       error('MEOrder: Invalid ''kind'' parameter!')
   end
end

