%  [D0, D1, logli] = MAPFromTrace(trace, orders, maxIter, stopCond, initial, result)
%  
%  Performs MAP fitting using the EM algorithm (ErCHMM, 
%  [1]_, [2]_).
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
%  initial : tuple of a vector and a matrix, shape(N,N), optional
%      The rate parameters of the Erlang distributions 
%      and the branch transition probability matrix to be
%      used initially. If not given, a default initial 
%      guess is determined and the algorithm starts from 
%      there.
%  result : {"vecmat", "matmat"}, optional
%      The result can be returned two ways. If "matmat" is
%      selected, the result is returned in the classical
%      representation of MAPs, thus the D0 and D1 matrices.
%      If "vecmat" is selected, the rate parameters of the
%      Erlang branches and the branch transition probability
%      matrix are returned. The default value is "matmat"
%  
%  Returns
%  -------
%  (D0, D1) : tuple of matrix, shape (M,M) and matrix, shape (M,M)
%      If the "matmat" result format is chosen, the function
%      returns the D0 and D1 matrices of the MAP
%  (lambda, P) : tuple of vector, length N and matrix, shape (M,M)
%      If the "vecmat" result format is chosen, the function
%      returns the vector of the Erlang rate parameters of 
%      the branches and the branch transition probability 
%      matrix
%  logli : double
%      The log-likelihood divided by the trace length
%      
%  Notes
%  -----
%  This procedure is quite slow in the supported 
%  mathematical frameworks. If the maximum speed is
%  needed, please use the multi-core optimized c++
%  implementation called SPEM-FIT_.
%  
%  .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit
%  
%  References
%  ----------
%  .. [1] Okamura, Hiroyuki, and Tadashi Dohi. Faster 
%         maximum likelihood estimation algorithms for 
%         Markovian arrival processes. Quantitative 
%         Evaluation of Systems, 2009. QEST'09. Sixth 
%         International Conference on the. IEEE, 2009.
%  
%  .. [2] Horváth, Gábor, and Hiroyuki Okamura. A Fast EM
%         Algorithm for Fitting Marked Markovian Arrival 
%         Processes with a New Special Structure. Computer
%         Performance Engineering. Springer Berlin 
%         Heidelberg, 2013. 119-133.

function [D0, D1, logli] = MAPFromTrace (trace, orders, maxIter, stopCond, initial, result)

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
        result = 'matmat';
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
    function printOrders (ord)
        for oi=1:length(ord)
            fprintf('%g',ord(oi));
            if oi<length(ord)
                fprintf(',');
            end
        end           
    end

    if length(orders)==1
        bestX = [];
        bestY = [];
        bestLogli = -inf;
        bestOrders= [];
        for br=2:orders
            allord = allorders (br, orders);
            for ox=1:length(allord)
                if BuToolsVerbose
                    fprintf('Trying orders ');
                    printOrders(allord{ox});
                    fprintf('...\n');
                end
                [oX,oY,l] = MAPFromTrace (trace, allord{ox}, maxIter, stopCond, initial, result);
                if l > bestLogli
                    bestX = oX;
                    bestY = oY;
                    bestLogli = l;
                    bestOrders = allord{ox};
                end
            end
        end
        D0 = bestX;
        D1 = bestY;
        logli = bestLogli;
        if BuToolsVerbose
            fprintf('Best solution: logli=%g, orders=', bestLogli);
            printOrders(bestOrders);
            fprintf('\n');
        end
        return;
    end
    
    function X = generatorFromErlangs (erll, erlo)
        X = zeros (sum(erlo));
        xx = 1;
        for ii=1:length(erll)
            X(xx:xx+erlo(ii)-1, xx:xx+erlo(ii)-1) = erll(ii)*(diag(ones(1,erlo(ii)-1),1)-diag(ones(1,erlo(ii))));
            xx = xx + erlo(ii);
        end
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
        P = ones(M,1)*alphav;
    elseif length(initial)==2
        lambda = initial{1};
        P = initial{2};
        alphav = DTMCSolve(P);
    else
        error('Invalid initial branch probability and rate vectors!');
    end

    Q = zeros(M, K);
    A = zeros(K, M);
    B = zeros(M, K);
    Ascale = zeros(1,K);
    Bscale = zeros(1,K);
    ologli = 1;
    logli = 0;
    steps = 1;
    t1 = clock();
    while abs((ologli-logli)/logli)>stopCond && steps<=maxIter
        ologli = logli;
        % E-step:
        % branch densities:
        for i=1:M
            Q(i,:) = ((lambda(i)*trace).^(orders(i)-1) / factorial(orders(i)-1) * lambda(i)) .* exp(-lambda(i)*trace);
        end
        % forward likelihood vectors:
        prev = alphav;
        scprev = 0;
        for k=1:K
            prev = prev*diag(Q(:,k))*P;
            scale = log2(sum(prev));
            prev = prev * 2^-scale;
            Ascale(k) = scprev + scale;
            A(k,:) = prev;
            scprev = Ascale(k);
        end
        Av = [alphav; A(1:end-1,:)];
        Ascalev = [0, Ascale(1:end-1)];
        % backward likelihood vectors:
        next = ones(M,1);
        scprev = 0;
        for k=K:-1:1
            next = diag(Q(:,k))*P*next;
            scale = log2(sum(next));
            next = next * 2^-scale;
            Bscale(k) = scprev + scale;
            B(:,k) = next;
            scprev = Bscale(k);
        end
        Bv = [B(:,2:end), ones(M,1)];
        Bscalev = [Bscale(2:end), 0];
        
        llh = alphav*B(:,1);
        logli = (log(llh) + Bscale(1) * log(2)) / K;
        illh = 1.0 / llh;

        % M-step:
        % Calculate new estimates for the parameters
        AB = Av.*B';
        nor = sum(AB,2);
        for m=1:M
            AB(:,m) = AB(:,m)./nor;
        end
        v1 = sum(AB,1);
        v2 = trace'*AB;
        alphav = v1/K;
        lambda = (orders.*v1 ./ v2)';
        
        Avv = Av.*Q';
        nor = illh*2.^(Ascalev+Bscalev-Bscale(1))';
        for m=1:M
            Avv(:,m) = Avv(:,m) .* nor;
        end       
        P = (Avv'*Bv').*P;
        for m=1:M
            P(m,:) = P(m,:) / sum(P(m,:));
        end        
        
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

    if strcmp(result,'vecmat')
        D0 = lambda;
        D1 = P;
    elseif strcmp(result,'matmat')
        D0 = generatorFromErlangs(lambda, orders);
        D1 = zeros(size(D0));
        indicesTo = [1 cumsum(orders(1:end-1))+1];
        indicesFrom = cumsum(orders);
        D1(indicesFrom,indicesTo) = diag(lambda)*P;
    end
end

