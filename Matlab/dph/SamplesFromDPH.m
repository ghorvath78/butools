%  x = SamplesFromDPH(alpha, A, K, prec)
%  
%  Generates random samples from a discrete phase-type 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the discrete phase-
%      type distribution.
%  A : matrix, shape (M,M)
%      The transition probability  matrix of the discrete phase-
%      type distribution.
%  K : integer
%      The number of samples to generate.
%  prec : double, optional
%      Numerical precision to check if the input phase-type
%      distribution is valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  x : vector, length(K)
%      The vector of random samples

function x = SamplesFromDPH(a,A,k)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDPHRepresentation(a,A)
        error('SamplesFromDPH: input isn''t a valid DPH representation!');
    end

    % auxilary variables
    N = length(a);
    cummInitial = cumsum(a);
    logp = log(diag(A));
    sojourn = 1./(1-diag(A));
    nextpr = diag(sojourn)*A;
    nextpr = nextpr - diag(diag(nextpr));
    nextpr = [nextpr, 1-sum(nextpr,2)];
    nextpr = cumsum(nextpr,2);
    
    x = zeros(k,1);
    for n=1:k
        time = 0;

        % draw initial distribution
        r = rand();
        state = 1;
        while cummInitial(state)<=r
            state=state+1;
        end

        % play state transitions
        while state<=N 
            time = time + 1 + floor(log(rand()) / logp(state));
            r = rand();
            nstate = 1;
            while nextpr(state,nstate)<=r
                nstate = nstate+1;
            end
            state = nstate;
        end
        x(n) = time;
    end
end
