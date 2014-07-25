%  [alpha, A, logli] = PHFromTrace(trace, orders, maxIter, stopCond, initial, result)
%  
%  Performs PH distribution fitting using the EM algorithm
%  (G-FIT, [1]_).
%  
%  Parameters
%  ----------
%  trace : column vector, length K
%      The samples of the trace
%  orders : list of int, length(N), or int
%      The length of the list determines the number of 
%      Erlang branches to use in the fitting method.
%      The entries of the list are the orders of the 
%      Erlang distributions. If this parameter is a 
%      single integer, all possible branch number - order
%      combinations are tested where the total number of 
%      states is "orders".
%  maxIter : int, optional
%      Maximum number of iterations. The default value is
%      200
%  stopCond : double, optional
%      The algorithm stops if the relative improvement of
%      the log likelihood falls below stopCond. The 
%      default value is 1e-7
%  initial : tuple of two vectors, optional
%      The initial values of the branch probabilities and
%      rate parameters is given by this tuple. If not 
%      given, a default initial guess is determined and 
%      the algorithm starts from there.
%  result : {"vecmat", "vecvec"}, optional
%      The result can be returned two ways. If "vecmat" is
%      selected, the result is returned in the classical
%      representation of phase-type distributions, thus the
%      initial vector and the generator matrix. 
%      If "vecvec" is selected, two vectors are returned, 
%      one holds the branch probabilities, and the second
%      holds the rate parameters of the Erlang branches.
%      The default value is "vecmat"
%  
%  Returns
%  -------
%  (alpha, A) : tuple of matrix, shape (1,M) and matrix, shape (M,M)
%      If the "vecmat" result format is chosen, the function
%      returns the initial probability vector and the
%      generator matrix of the phase type distribution.
%  (pi, lambda) : tuple of vector, length N and vector, length N
%      If the "vecvec" result format is chosen, the function
%      returns the vector of branch probabilities and the
%      vector of branch rates in a tuple.
%  logli : double
%      The log-likelihood divided by the trace length
%      
%  Notes
%  -----
%  This procedure is quite fast in the supported 
%  mathematical frameworks. If the maximum speed is
%  needed, please use the multi-core optimized c++
%  implementation called SPEM-FIT_.
%  
%  .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit
%  
%  References
%  ----------
%  .. [1] Thummler, Axel, Peter Buchholz, and Miklós Telek.
%         A novel approach for fitting probability 
%         distributions to real trace data with the EM 
%         algorithm. Dependable Systems and Networks, 2005.

function [alpha, A, logli] = PHFromTrace (trace, orders, maxIter, stopCond, initial, result)

    if ~exist('maxIter','var') || isempty(maxIter)
        maxIter = 200;
    end
    if ~exist('stopCond','var') || isempty(stopCond)
        stopCond = 1e-7;
    end
    if ~exist('initial','var')
        initial = {};
    end
    if ~exist('result','var') || isempty(result)
        result = 'vecmat';
    end    
    global BuToolsVerbose;

    function o = allorders (branches, sumorders)
        if branches==1
            o={sumorders};
        else
            o = {};
            for ii=1:sumorders-branches+1
                x = allorders (branches-1, sumorders-ii);
                for j=1:length(x)
                    xt = sort([x{j} ii]);
                    % check if we have it already
                    found=false;
                    for ok=1:length(o)
                        if o{ok}==xt
                            found=true;
                            break;
                        end
                    end
                    if ~found
                        o{length(o)+1}=xt;
                    end
                end
            end
        end
    end

    if length(orders)==1
        bestAlpha = [];
        bestA = [];
        bestLogli = -inf;
        bestOrders= [];
        for br=2:orders
            allord = allorders (br, orders);
            for ox=1:length(allord)
                [a,A,l] = PHFromTrace (trace, allord{ox}, maxIter, stopCond, initial, result);
                if l > bestLogli
                    bestAlpha = a;
                    bestA = A;
                    bestLogli = l;
                    bestOrders = allord{ox};
                end
            end
        end
        alpha = bestAlpha;
        A = bestA;
        logli = bestLogli;
        if BuToolsVerbose
            fprintf('Best solution: logli=%g, orders=', bestLogli);
            for i=1:length(bestOrders)
                fprintf('%g',bestOrders(i));
                if i<length(bestOrders)
                    fprintf(',');
                else
                    fprintf('\n');
                end
            end           
        end
        return;
    end
    
    M = length(orders);
    K = length(trace);

    % initial alpha and lambda is such that the mean is matched
    if isempty(initial)
        alphav = ones(1,M) / M;
        lambda = diag(orders) * (1:M)';
        trm = sum(trace)/length(trace);
        inim = sum(alphav ./ (1:M));
        lambda = lambda * inim / trm;
    elseif length(initial)==2
        if length(initial{1})==M && length(initial{2})==M
            alphav = reshape(initial{1},1,M);
            lambda = reshape(initial{2},M,1);
        else
            error('The length of the initial branch probability and rate vectors is not consistent with the length of the orders vector!');
        end
    else
        error('Invalid initial branch probability and rate vectors!');
    end

    Q = zeros(M, K);
    ologli = 1;
    logli = 0;
    steps = 1;
    t1 = clock();
    while abs((ologli-logli)/logli)>stopCond && steps<=maxIter
        ologli = logli;
        % E-step:
       for i=1:M
           Q(i,:) = (alphav(i)*(lambda(i)*trace).^(orders(i)-1) / factorial(orders(i)-1) * lambda(i)) .* exp(-lambda(i)*trace);
       end
       logli = sum(log(sum(Q,1))) / K;
       nor = sum(Q,1);
       for i=1:M
           Q(i,:) = Q(i,:) ./ nor;
       end
       % M-step:
       v1 = sum(Q,2);
       v2 = Q*trace;
       alphav = v1'/K;
       lambda = diag(orders)*v1 ./ v2;
       steps = steps+1;
       if BuToolsVerbose && etime(clock(),t1)>2
            fprintf('Num of iterations: %d, logli: %g\n', steps, logli);
           t1 = clock();
       end
    end

    if BuToolsVerbose
        fprintf('Num of iterations: %d, logli: %g\n', steps, logli);
        fprintf('EM algorithm terminated. (orders=');
        for i=1:M
            fprintf('%g',orders(i));
            if i<M
                fprintf(',');
            else
                fprintf(')\n');
            end
        end
    end

    if strcmp(result,'vecvec')
        alpha = alphav;
        A = lambda;
    elseif strcmp(result,'vecmat')
        N = sum(orders);
        alpha = zeros(1,N);
        A = zeros (N);
        ix = 1;
        for i=1:M
            alpha(ix) = alphav(i);
            A(ix:ix+orders(i)-1, ix:ix+orders(i)-1) = lambda(i)*(diag(ones(1,orders(i)-1),1)-diag(ones(1,orders(i))));
            ix = ix + orders(i);
        end
    end
end

