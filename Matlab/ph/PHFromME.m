%  [beta, B] = PHFromME(alpha, A, precision)
%  
%  Obtains a Markovian representation of a matrix 
%  exponential distribution of the same size, if possible.
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-exponential
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      distribution.
%  precision : double, optional
%      A representation is considered to be a Markovian one
%      if it is closer than the precision. The default value
%      is 1e-14.
%  
%  Returns
%  -------
%  beta : vector, shape (1,M)
%      The initial probability vector of the Markovian 
%      monocyclic representation
%  B : matrix, shape (M,M)
%      Transient generator matrix of the Markovian 
%      monocyclic representation
%  
%  References
%  ----------
%  .. [1] G Horváth, M Telek, "A minimal representation of 
%         Markov arrival processes and a moments matching 
%         method," Performance Evaluation 64:(9-12) pp. 
%         1153-1168. (2007)

function [beta, B] = PHFromME (alpha, A, precision)

    function nrep = transfun (orep, B)
        nrep = {orep{1}*B, inv(B)*orep{2}*B};
    end
        
    function d = evalfun (orep, k)
        if nargin<2
            k = 0;
        end
        ao = orep{1};
        Ao = orep{2};
        av = sum(-Ao,2);
        Ad = Ao - diag(diag(Ao));
        if rem(k,2) == 0
            d = -min([min(ao), min(av), min(min(Ad))]);
        else
            d = -sum(ao(ao<0)) - sum(av(av<0)) - sum(sum(Ad(Ad<0)));
        end
    end

    if ~exist('precision','var')
        precision = 1e-14;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMERepresentation(alpha, A)
        error('PHFromME: Input isn''t a valid ME distribution!');
    end

    nrep = FindMarkovianRepresentation ({alpha, A}, @transfun, @evalfun, precision);
    beta = nrep{1};
    B = nrep{2};
end

