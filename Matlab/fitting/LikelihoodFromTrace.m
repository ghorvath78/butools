%  logli = LikelihoodFromTrace(trace, X, Y, prec)
%  
%  Evaluates the log-likelihood of a trace with the given PH
%  distribution or MAP. The result is divided by the length
%  of the trace.
%  
%  If X is a row vector, than (X,Y) is interpreted as a PH
%  distribution, otherwise (X,Y) is considered to be a MAP.
%  
%  Parameters
%  ----------
%  trace : column vector, length K
%      The samples of the trace
%  X : matrix, shape (1,M) or (M,M)
%      If X is a row vector, it is the initial probability
%      vector of the PH distribution. If X is a square
%      matrix, it is interpreted as the D0 matrix of a MAP
%  Y : matrix, (M,M)
%      If X is a row vector, Y is the transient generator
%      of the PH distribution. If X is a square matrix, Y
%      is interpreted as the D1 matrix of a MAP
%  prec : double, optional
%      Numerical precision used by the randomization. The
%      default value is 1e-14.
%  
%  Returns
%  -------
%  logli : double
%      The log likelihood divided by the size of the trace
%      
%  Notes
%  -----
%  The procedure is much faster with PH distributions.

function logli = LikelihoodFromTrace (trace, X, Y, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if size(X,1)==1

        % We have a PH distribution. We can sort it to make the computation
        % faster
        alpha = X;
        A = Y;
        trace = sort(trace);
        lambda = max(abs(diag(A)));
        P = A/lambda + eye(size(A));
        a = sum(-A,2);
        eps = max(prec, 10^(log10(prec) + log10(lambda)));
        lpoi = -lambda*trace;
        trace = log(trace);
        poi =  exp(lpoi);
        spoi = poi;
        fx = poi*(alpha*a);
        k = 1;
        first = 1;
        coeffv = alpha;
        maxIter = 10000;
        while ~isempty(first) && k<maxIter
            coeffv = coeffv * P;
            lpoi(first:end) = lpoi(first:end) + log(lambda) + trace(first:end) - log(k);
            poi(first:end) = exp(lpoi(first:end));
            spoi(first:end) = spoi(first:end) + poi(first:end);
            fx(first:end) = fx(first:end) + poi(first:end) * (coeffv * a);
            k = k + 1;
            first = first + find(spoi(first:end)<1-eps,1,'first') - 1;
        end
        logli = sum(log(fx))/length(trace);
    else        
        D0 = X;
        D1 = Y;
        N = size(D0,1);
        L = length(trace);
        
        % first we calculate matrix e^(D0*x(i))*D1 for each sample
        [trace,ix] = sort(reshape(trace,[],1));
        lambda = max(abs(diag(D0)));
        P = D0/lambda + eye(size(D0));
        eps = max(prec, 10^(log10(prec) + log10(lambda)));
        lpoi = -lambda*trace;
        trace = log(trace);
        poi =  exp(lpoi);
        spoi = poi;
        coeffv = D1;
        fx = kron(poi,coeffv);
        k = 1;
        first = 1;
        maxIter = 10000;
        while ~isempty(first) && k<maxIter
            coeffv = P * coeffv;
            lpoi(first:end) = lpoi(first:end) + log(lambda) + trace(first:end) - log(k);
            poi(first:end) = exp(lpoi(first:end));
            spoi(first:end) = spoi(first:end) + poi(first:end);
            fx((first-1)*N+1:end,:) = fx((first-1)*N+1:end,:) + kron(poi(first:end),coeffv);
            k = k + 1;
            first = first + find(spoi(first:end)<1-eps,1,'first') - 1;
        end
        alpha = DTMCSolve (inv(-D0)*D1);
        l = alpha;
        sc = 0;
        [~,ixrev]=sort(ix);
        for i=1:L
            l = l*fx((ixrev(i)-1)*N+1:ixrev(i)*N,:);
            if rem(i,10)==0
                % sometimes we need to rescale the results to avoid "nan"s
                scale = ceil(log2(sum(l)));            
                if scale>1
                    l = l/2^scale;
                    sc = sc+scale;
                end
                if scale<-10
                    scale = scale+10;
                    l = l/2^scale;
                    sc = sc+scale;
                end
            end
        end
        logli = (log(sum(l))+sc*log(2)) / length(trace);
    end
end