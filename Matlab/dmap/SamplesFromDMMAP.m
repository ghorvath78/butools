%  x = SamplesFromDMMAP(D, K, prec)
%  
%  Generates random samples from a discrete marked 
%  Markovian arrival process.
%  
%  Parameters
%  ----------
%  D : list of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the DMMAP
%  K : integer
%      The number of samples to generate.
%  prec : double, optional
%      Numerical precision to check if the input DMMAP is
%      valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  x : matrix, shape(K,2)
%      The random samples. Each row consists of two 
%      columns: the (discrete) inter-arrival time and the
%      type of the arrival.        

function x = SamplesFromDMMAP(D,k,initial)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMMAPRepresentation(D)
        error('SamplesFromDMMAP: input isn''t a valid DMMAP representation!');
    end

    N = size(D{1},1);
    
    if ~exist('initial','var') || isempty(initial)
        % draw initial state according to the stationary distribution
        stst = MarginalDistributionFromDMMAP(D);
        cummInitial = cumsum(stst);
        r = rand();
        state = 1;
        while cummInitial(state)<=r
            state=state+1;
        end
    else
        state = initial;
    end

    % auxilary variables
    sojourn = 1./(1-diag(D{1}));
    logp = log(diag(D{1}));
    nextpr = diag(sojourn)*D{1};
    nextpr = nextpr - diag(diag(nextpr));
    for i=2:length(D)
        nextpr=[nextpr,diag(sojourn)*D{i}];
    end
    nextpr = cumsum(nextpr,2);
    
    if length(D)>2
        x = zeros(k,2);
    else
        x = zeros(k,1);
    end        
    for n=1:k
        time = 0;

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
        if length(D)>2
            x(n,1) = time;
            x(n,2) = ceil(state/N)-1;
        else
            x(n) = time;
        end
        state = rem(state-1,N)+1;
    end
end
