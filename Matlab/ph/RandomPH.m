%  [alpha, A] = RandomPH(order, mean, zeroEntries, maxTrials, prec)
%  
%  Returns a random phase-type distribution with a given 
%  order.
%  
%  Parameters
%  ----------
%  order : int
%      The size of the phase-type distribution
%  mean : double, optional
%      The mean of the phase-type distribution 
%  zeroEntries : int, optional
%      The number of zero entries in the initial vector, 
%      generator matrix and closing vector
%  maxTrials : int, optional
%      The maximum number of trials to find a proper PH 
%      (that has an irreducible phase process and none of 
%      its parameters is all-zero). The default value is 
%      1000.
%  prec : double, optional
%      Numerical precision for checking the irreducibility.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  alpha : vector, shape (1,M)
%      The initial probability vector of the phase-type 
%      distribution.
%  A : matrix, shape (M,M)
%      The transient generator matrix of the phase-type 
%      distribution.
%  
%  Notes
%  -----
%  If the procedure fails, try to increase the 'maxTrials'
%  parameter.   

function [alpha, A] = RandomPH(order, mean, zeroEntries, maxTrials, prec)

    if ~exist('zeroEntries','var')
        zeroEntries = 0;
    end

    if ~exist('prec','var')
        prec = 1e-7;
    end

    if ~exist('mean','var')
        mean = 1;
    end

    if ~exist('maxTrials','var')
        maxTrials = 1000;
    end

    if zeroEntries > (order+1)*(order-1)
        error('RandomPH: Too many zero entries requested! Try to decrease the zero entries number!');
    end

    % distribute the zero entries among the rows
    function o = allZeroDistr (states, zeros)
        if states==1
            o = zeros;
        else
            o = [];
            for iz=0:zeros
                x = allZeroDistr (states-1, zeros-iz);
                for jz=1:size(x,1)
                    xt = sort([x(jz,:), iz]);
                    % check if we have it already
                    found = 0;
                    for kz=1:size(o,1)
                        if sum((o(kz,:)-xt).^2)==0
                            found = 1;
                            break;
                        end
                    end
                    if ~found
                        o = [o; xt];
                    end
                end
            end
        end
    end
        
    zeroDistr = allZeroDistr(order, zeroEntries);   

    trials = 1;
    while trials<maxTrials
        % select a configuration from zeroDistr: it is a list describing the zero entries in each row
        zdix = randperm(size(zeroDistr,1));
        for k=1:size(zeroDistr,1)
            zDistr = zeroDistr(zdix(k),:);           
            B = zeros(order,order+2);
            for i=1:order
                rp = randperm(order+1);
                a = zeros(1,order+1);
                for j=1:order+1-zDistr(i)
                    a(rp(j)) = rand();
                end
                B(i,1:i-1) = a(1:i-1);
                B(i,i+1:end) = a(i:end);
            end
            % construct PH parameters
            A = B(:,1:order);
            a = B(:,order+2);
            A = A - diag(sum(A,2)+a);
            alpha = B(:,order+1)';
            % check if it is a proper PH (irreducible phase process & no full zero matrix)
            if all(all(A==0.0)) || all(alpha==0.0) || all(a==0.0)
                continue
            end
            alpha = alpha / sum(alpha);
            D = A + a*alpha;
            if rank(D) == order-1
                pi = CTMCSolve(D);
                if min(abs(pi)) > prec
                    % scale to the mean value
                    m = MomentsFromPH (alpha, A, 1);
                    A = A * m(1) / mean;
                    return;
                end
            end
            trials = trials + 1;
        end
    end
    error ('No feasible random PH found with such many zero entries!');
end
