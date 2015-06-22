%  mrep = FindMarkovianRepresentation(rep, @transfun, @evalfunc, precision)
%  
%  Obtains a Markovian representation from a non-Markovian 
%  one while keeping the size the same, by applying a series 
%  of elementary transformations.
%  
%  Parameters
%  ----------
%  rep : tuple of matrices
%      The initial non-Markovian representation
%      (initial vector and generator of a PH, matrices of a
%      MAP, or a MMAP, etc.)
%  transfun : callable
%      A function that transforms the representation using 
%      the given similarity transformation matrix
%  evalfunc : callable
%      A function that returns how far the representation is
%      from the Markovian one
%  precision : double
%      A representation is considered to be a Markovian one
%      if it is closer than the precision. The default value
%      is 1e-7
%      
%  Returns
%  -------
%  mrep : tuple of matrices
%      The Markovian representation, if found. If not found,
%      the closest one is returned.
%  
%  Notes
%  -----
%  This function should not be called directly.
%  It is used by 'PHFromME', 'MAPFromRAP', etc. functions.
%  
%  References
%  ----------
%  .. [1]  G Horváth, M Telek, "A minimal representation of 
%          Markov arrival processes and a moments matching 
%          method," Performance Evaluation 64:(9-12) 
%          pp. 1153-1168. (2007)

function nrep = FindMarkovianRepresentation (rep, transfun, evalfun, precision)

    if ~exist('precision','var')
        precision = 1e-7;
    end

    function [bestrep, bestdist] = elementary (erep, b, k, evalfun, transfun)
        bestdist = evalfun (erep, k);
        bestrep = erep;
        repSize = size(erep{1},2);
        for i=1:repSize
            for j=1:repSize
                if i~=j
                    % create elementary transformation matrix with +b
                    B = eye(repSize);
                    B(i,j) = b;
                    B(i,i) = 1.0 - b;
                    % apply similarity transform
                    newrep = transfun (erep, B);
                    newdist = evalfun (newrep, k);
                    % store result if better
                    if newdist < bestdist
                        bestrep = newrep;
                        bestdist = newdist;
                    end
                    % create elementary transformation matrix with -b
                    B = eye(repSize);
                    B(i,j) = -b;
                    B(i,i) = 1.0 + b;
                    % apply similarity transform
                    newrep = transfun (erep, B);
                    newdist = evalfun (newrep, k);
                    % store result if better
                    if newdist < bestdist
                        bestrep = newrep;
                        bestdist = newdist;
                    end
                end
            end
        end
    end

    function [bestrep, lastdist] = minimize (orep, iters, b, k, evalfun, transfun)
        lastdist = evalfun (orep, k);
        bestrep = orep;
        for i=1:iters
            [orep, dist] = elementary (orep, b, k, evalfun, transfun);
            if dist >= lastdist
                break;
            else
                lastdist = dist;
                bestrep = orep;
            end
        end
    end

    if evalfun(rep) < precision
        nrep = rep;
        return;
    end
    
    nrep = rep;
    M = size(nrep{1},2);
    b = 0.5;
    odist = inf;
    while b>precision/2
        for m=1:M*M
            for k=1:4
                [nrep, ddist] = minimize (nrep, M*M, b, k, evalfun, transfun);
                if ddist < precision
                    return;
                end
            end
            if odist <= ddist
                break;
            end
            odist = ddist;
        end
        b = b / 2.0;
    end
end
