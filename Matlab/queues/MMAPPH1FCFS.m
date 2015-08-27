%  Ret = MMAPPH1FCFS(D, sigma, S, ...)
%  
%  Returns various performane measures of a MMAP[K]/PH[K]/1 
%  first-come-first-serve queue, see [1]_.
%  
%  Parameters
%  ----------
%  D : list of matrices of shape (N,N), length (K+1)
%      The D0...DK matrices of the arrival process.
%  sigma : list of row vectors, length (K)
%      The list containing the initial probability vectors of the service
%      time distributions of the various customer types. The length of the
%     vectors does not have to be the same.
%  S : list of square matrices, length (K)
%      The transient generators of the phase type distributions representing
%      the service time of the jobs belonging to various types.
%  further parameters : 
%      The rest of the function parameters specify the options
%      and the performance measures to be computed.
%  
%      The supported performance measures and options in this 
%      function are:
%  
%      +----------------+--------------------+----------------------------------------+
%      | Parameter name | Input parameters   | Output                                 |
%      +================+====================+========================================+
%      | "ncMoms"       | Number of moments  | The moments of the number of customers |
%      +----------------+--------------------+----------------------------------------+
%      | "ncDistr"      | Upper limit K      | The distribution of the number of      |
%      |                |                    | customers from level 0 to level K-1    |
%      +----------------+--------------------+----------------------------------------+
%      | "stMoms"       | Number of moments  | The sojourn time moments               |
%      +----------------+--------------------+----------------------------------------+
%      | "stDistr"      | A vector of points | The sojourn time distribution at the   |
%      |                |                    | requested points (cummulative, cdf)    |
%      +----------------+--------------------+----------------------------------------+
%      | "stDistrME"    | None               | The vector-matrix parameters of the    |
%      |                |                    | matrix-exponentially distributed       |
%      |                |                    | sojourn time distribution              |
%      +----------------+--------------------+----------------------------------------+
%      | "stDistrPH"    | None               | The vector-matrix parameters of the    |
%      |                |                    | matrix-exponentially distributed       |
%      |                |                    | sojourn time distribution, converted   |
%      |                |                    | to a continuous PH representation      |
%      +----------------+--------------------+----------------------------------------+
%      | "prec"         | The precision      | Numerical precision used as a stopping |
%      |                |                    | condition when solving the Riccati     |
%      |                |                    | equation                               |
%      +----------------+--------------------+----------------------------------------+
%      | "classes"      | Vector of integers | Only the performance measures          |
%      |                |                    | belonging to these classes are         |
%      |                |                    | returned. If not given, all classes    |
%      |                |                    | are analyzed.                          |
%      +----------------+--------------------+----------------------------------------+
%      
%      (The quantities related to the number of customers in 
%      the system include the customer in the server, and the 
%      sojourn time related quantities include the service 
%      times as well)
%  
%  Returns
%  -------
%  Ret : list of the performance measures
%      Each entry of the list corresponds to a performance 
%      measure requested. Each entry is a matrix, where the
%      columns belong to the various job types.
%      If there is just a single item, 
%      then it is not put into a list.
%  
%  References
%  ----------
%  .. [1] Qiming He, "Analysis of a continuous time 
%         SM[K]/PH[K]/1/FCFS queue: Age process, sojourn times,
%         and queue lengths", Journal of Systems Science and 
%         Complexity, 25(1), pp 133-155, 2012.

function varargout = MMAPPH1FCFS(D, sigma, S, varargin)
    
    K = length(D)-1;

    % parse options
    precision = 1e-14;
    classes = 1:K;
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'prec')
            precision = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif strcmp(varargin{i},'classes')
            classes = varargin{i+1};
            eaten = [eaten, i, i+1];
        end
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMMAPRepresentation(D)
        error('MMAPPH1FCFS: The arrival process is not a valid MMAP representation!');
    end
    
    if BuToolsCheckInput
        for k=1:K
            if ~CheckPHRepresentation(sigma{k},S{k})
                error('MMAPPH1FCFS: the vector and matrix describing the service times is not a valid PH representation!');
            end
        end
    end

    % Various properties of the arrival and service processes
    D0 = D{1};
    N = size(D0,1);
    Ia = eye(N);
    Da = zeros(N);
    for q=1:K
        Da = Da + D{q+1};
    end
    theta = CTMCSolve(D0+Da);
    beta = cell(1,K);
    lambda = zeros(K,1);
    mu = zeros(K,1);
    Nsk = zeros(1,K);
    ro = 0;
    for k=1:K
        lambda(k) = sum(theta*D{k+1});
        beta{k} = CTMCSolve(S{k} - sum(S{k},2)*sigma{k});
        mu(k) = sum(beta{k}*(-S{k}));
        Nsk(k) = size(S{k},1);
        ro = ro + lambda(k)/mu(k);
    end
    alpha = theta*Da/sum(lambda);
    D0i = inv(-D0);

    Sa = S{1};
    sa = cell(1,K);
    ba = cell(1,K);
    sa{1} = sigma{1};
    ba{1} = beta{1};
    sv{1} = -sum(S{1},2);
    Pk = cell(1,K);
    Pk{1} = D0i*D{2};
    for q=2:K
        sa{q} = zeros(1,length(sigma{1}));
        ba{q} = zeros(1,length(beta{1}));
        sv{q} = zeros(length(sigma{1}),1);
        Pk{q} = D0i*D{q+1};
    end
    for k=2:K
        Sa = blkdiag (Sa, S{k});
        for q=1:K
            if q==k
                sa{q} = [sa{q} sigma{k}];
                ba{q} = [ba{q} beta{k}];
                sv{q} = [sv{q}; -sum(S{k},2)];
            else
                sa{q} = [sa{q} zeros(size(sigma{k}))];
                ba{q} = [ba{q} zeros(size(beta{k}))];
                sv{q} = [sv{q}; zeros(length(sigma{k}),1)];
            end
        end
    end
    P = D0i*Da;
    iVec = kron(D{2},sa{1});
    for k=2:K
        iVec = iVec + kron(D{k+1},sa{k});
    end
    Ns = size(Sa,1);
    Is = eye(Ns);
    
    % step 1. solve the age process of the queue
    % ==========================================

    % solve Y0 and calculate T
    Y0 = FluidFundamentalMatrices (kron(Ia,Sa), kron(Ia,-sum(Sa,2)), iVec, D0, 'P', precision);
    T = kron(Ia,Sa) + Y0 * iVec;

    % calculate pi0 and v0
    pi0 = zeros(1,size(T,1));
    for k=1:K
        pi0 = pi0 + kron(theta*D{k+1},ba{k}/mu(k));
    end
    pi0 = - pi0 * T;

    iT = inv(-T);
    oa = ones(N,1);

    % step 2. calculate performance measures
    % ======================================
    Ret = {};
    for k=classes
        argIx = 1;
        clo = iT*kron(oa,sv{k});
        while argIx<=length(varargin)
            if any(ismember(eaten, argIx))
                argIx = argIx + 1;
                continue;
            elseif strcmp(varargin{argIx},'stMoms') 
                numOfSTMoms = varargin{argIx+1};
                rtMoms = zeros(1,numOfSTMoms);
                for m=1:numOfSTMoms
                    rtMoms(m) = factorial(m) * pi0 * iT^m * clo / (pi0*clo);
                end
                Ret{end+1} = rtMoms;
                argIx = argIx + 1;
            elseif strcmp(varargin{argIx},'stDistr') 
                stCdfPoints = varargin{argIx+1};
                cdf = [];
                for t=stCdfPoints
                    pr = 1 - pi0 * expm(T*t) * clo / (pi0*clo);
                    cdf = [cdf, pr];
                end
                Ret{end+1} = cdf;
                argIx = argIx + 1;
            elseif strcmp(varargin{argIx},'stDistrME')
                Bm = SimilarityMatrixForVectors(clo/(pi0*clo),ones(N*Ns,1));
                Bmi = inv(Bm);
                A = Bm * T * Bmi;
                alpha = pi0 * Bmi;
                Ret{end+1} = alpha;
                Ret{end+1} = A;
            elseif strcmp(varargin{argIx},'stDistrPH')
                vv = pi0*iT;
                ix = 1:N*Ns;
                nz = ix(vv>precision);
                delta = diag(vv(nz));   
                cl = -T*clo/(pi0*clo);
                alpha = cl(nz,:)'*delta;
                A = inv(delta)*T(nz,nz)'*delta;
                Ret{end+1} = alpha;
                Ret{end+1} = A;
            elseif strcmp(varargin{argIx},'ncDistr')
                numOfQLProbs = varargin{argIx+1};
                argIx = argIx + 1;
                values = zeros(1,numOfQLProbs);
                jm = zeros(Ns,1);
                jm(sum(Nsk(1:k-1))+1:sum(Nsk(1:k)),:) = 1;
                jmc = ones(Ns,1);
                jmc(sum(Nsk(1:k-1))+1:sum(Nsk(1:k)),:) = 0;
                LmCurr = lyap(T, kron(D0+Da-D{k+1},Is), eye(N*Ns));
                values(1) = 1-ro+pi0*LmCurr*kron(oa,jmc);
                for i=1:numOfQLProbs-1
                    LmPrev = LmCurr;
                    LmCurr = lyap(T, kron(D0+Da-D{k+1},Is), LmPrev*kron(D{k+1},Is));
                    values(i+1) = pi0*LmCurr*kron(oa,jmc) + pi0*LmPrev*kron(oa,jm);
                end                
                Ret{end+1} = values;
            elseif strcmp(varargin{argIx},'ncMoms')
                numOfQLMoms = varargin{argIx+1};
                jm = zeros(Ns,1);
                jm(sum(Nsk(1:k-1))+1:sum(Nsk(1:k)),:) = 1;
                ELn = {lyap(T, kron(D0+Da,Is), eye(N*Ns))};
                qlMoms = zeros(1,numOfQLMoms);                
                for n=1:numOfQLMoms
                    bino = 1;
                    Btag = zeros(N*Ns);
                    for i=0:n-1
                        Btag = Btag + bino * ELn{i+1};
                        bino = bino * (n-i) / (i+1);
                    end
                    ELn{n+1} = lyap(T, kron(D0+Da,Is), Btag*kron(D{k+1},Is));
                    qlMoms(n) = sum(pi0*ELn{n+1}) + pi0*Btag*kron(oa,jm);
                end
                Ret{end+1} = qlMoms;
                argIx = argIx + 1;
            else
                error (['MMAPPH1FCFS: Unknown parameter ' varargin{argIx}])
            end
            argIx = argIx + 1;
        end
    end   
    varargout = Ret;
end