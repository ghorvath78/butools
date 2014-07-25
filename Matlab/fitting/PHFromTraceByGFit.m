function [alpha, A] = PHFromTraceByGFit (trace, orders)

    M = length(orders);
    K = length(trace);

    % initial alpha and lambda is such that the mean is matched
    alphav = ones(1,M) / M;
    lambda = diag(orders) * (1:M)';
    trm = sum(trace)/length(trace);
    inim = sum(alphav ./ (1:M));
    lambda = lambda * inim / trm;

    Q = zeros(M, K);
    ologli = 1;
    logli = 0;
    while abs(ologli-logli)>1e-7
        ologli = logli;
        % E-step:
       for i=1:M
           Q(i,:) = (alphav(i)*(lambda(i)*trace).^(orders(i)-1) / factorial(orders(i)-1) * lambda(i)) .* exp(-lambda(i)*trace);
       end
       logli = sum(log(sum(Q,1))) / K
       nor = sum(Q,1);
       for i=1:M
           Q(i,:) = Q(i,:) ./ nor;
       end
       % M-step:
       v1 = sum(Q,2);
       v2 = Q*trace;
       alphav = v1'/K;
       lambda = diag(orders)*v1 ./ v2;
    end

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

