function varargout = MMAPPH1PRPR(D, sigma, S, varargin)
    
    K = length(D)-1;

    % parse options
    erlMaxOrder = 200;
    precision = 1e-14;
    classes = 1:K;
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'erlMaxOrder')
            erlMaxOrder = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif strcmp(varargin{i},'prec')
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

    if BuToolsCheckInput && ~CheckMMAPRepresentation(D,precision)
        error('MMAPPH1PRPR: The arrival process is not a valid MMAP representation!');
    end
    
    if BuToolsCheckInput
        for k=1:K
            if ~CheckPHRepresentation(sigma{k},S{k},precision)
                error('MMAPPH1PRPR: the vector and matrix describing the service times is not a valid PH representation!');
            end
        end
    end

    % some preparation
    D0 = D{1};
    N = size(D0,1);
    I = eye(N);
    sD = zeros(N);
    for i=1:K+1
        sD = sD + D{i};
    end      
    
    s = cell(length(S));
    M = zeros(1,K);
    for i=1:K
        s{i} = sum(-S{i},2);
        M(i) = M(i) + length(sigma{i});
    end   
    
    Ret = {};
    for k=classes
       
        % step 1. solution of the workload process of the system
        % ======================================================
        sM = sum(M(k:K));
        Qwmm = D0;
        for i=1:k-1
            Qwmm = Qwmm + D{i+1};
        end
        Qwpm = zeros(N*sM, N);
        Qwmp = zeros(N, N*sM);
        Qwpp = zeros(N*sM);    
        kix = 1;
        for i=k:K
            Qwmp(:,kix:kix+N*M(i)-1) = kron(D{i+1}, sigma{i});
            Qwpm(kix:kix+N*M(i)-1,:) = kron(I,s{i});
            Qwpp(kix:kix+N*M(i)-1,kix:kix+N*M(i)-1) = kron(I,S{i});
            kix = kix + N*M(i);
        end

        % calculate fundamental matrices
        [Psiw, Kw, Uw] = FluidFundamentalMatrices (Qwpp, Qwpm, Qwmp, Qwmm, 'PKU', precision);
        
        % calculate boundary vector
        Ua = ones(N,1) + 2*sum(Qwmp*inv(-Kw),2);   
        pm = linsolve ([Uw,Ua]', [zeros(1,N),1]')';
      
        Bw = zeros(N*sM, N);
        Bw(1:N*M(k),:) = kron(I,s{k});
        kappa = pm*Qwmp / sum(pm*Qwmp*inv(-Kw)*Bw);
        
        if k<K            
            % step 2. construct fluid model for the remaining sojourn time process
            % ====================================================================
            % (for each class except the highest priority)
            Qsmm = D0;
            for i=1:k
                Qsmm = Qsmm + D{i+1};
            end
            Np = size(Kw,1);
            Qspm = zeros(Np+N*sum(M(1:k-1)), N);
            Qsmp = zeros(N, Np+N*sum(M(1:k-1)));
            Qspp = zeros(Np+N*sum(M(1:k-1)));
            Qspp(1:Np,1:Np) = Kw;
            Qspm(1:Np,1:N) = Bw;
            kix = Np+1;
            for i=k+1:K
                Qsmp(:,kix:kix+N*M(i)-1) = kron(D{i+1}, sigma{i});
                Qspm(kix:kix+N*M(i)-1,:) = kron(I,s{i});
                Qspp(kix:kix+N*M(i)-1,kix:kix+N*M(i)-1) = kron(I,S{i});
                kix = kix + N*M(i);
            end           
            inis = [kappa, zeros(1,N*sum(M(k+1:K)))];
            Psis = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, Qsmm, 'P', precision);            
            
            % step 3. calculate the performance measures
            % ==========================================   
            retIx = 1;
            argIx = 1;
            while argIx<=length(varargin)
                res = [];
                if any(ismember(eaten, argIx))
                    argIx = argIx + 1;
                    continue;
                elseif strcmp(varargin{argIx},'stMoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = varargin{argIx+1};
                    Pn = {Psis};
                    rtMoms = zeros(numOfSTMoms,1);
                    for n=1:numOfSTMoms
                        A = Qspp + Psis*Qsmp;
                        B = Qsmm + Qsmp*Psis;
                        C = -2*n*Pn{n};
                        bino = 1;
                        for i=1:n-1
                            bino = bino * (n-i+1) / i;
                            C = C + bino * Pn{i+1}*Qsmp*Pn{n-i+1};
                        end
                        P = lyap(A, B, C);
                        Pn{n+1} = P;
                        rtMoms(n) = sum(inis*P*(-1)^n) / 2^n;
                    end
                    res = rtMoms;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'stDistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = varargin{argIx+1};
                    for t=stCdfPoints
                        L = erlMaxOrder;
                        lambda = L/t/2;
                        Psie = FluidFundamentalMatrices (Qspp-lambda*eye(size(Qspp)), Qspm, Qsmp, Qsmm-lambda*eye(size(Qsmm)), 'P', precision);
                        Pn = {Psie};
                        pr = sum(inis*Psie);
                        for n=1:L-1
                            A = Qspp + Psie*Qsmp - lambda*eye(size(Qspp));
                            B = Qsmm + Qsmp*Psie - lambda*eye(size(Qsmm));
                            C = 2*lambda*Pn{n};
                            for i=1:n-1
                                C = C + Pn{i+1}*Qsmp*Pn{n-i+1};
                            end
                            P = lyap(A, B, C);
                            Pn{n+1} = P;
                            pr = pr + sum(inis*P);
                        end
                        res = [res; pr];
                    end
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'qlMoms')
                    % MOMENTS OF THE NUMBER OF JOBS
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLMoms = varargin{argIx+1};
                    % first calculate it at departure instants
                    QLDPn = {Psis};
                    dqlMoms = zeros(numOfQLMoms,1);
                    for n=1:numOfQLMoms
                        A = Qspp + Psis*Qsmp;
                        B = Qsmm + Qsmp*Psis;
                        C = n*QLDPn{n}*D{k+1};
                        bino = 1;
                        for i=1:n-1
                            bino = bino * (n-i+1) / i;
                            C = C + bino * QLDPn{i+1}*Qsmp*QLDPn{n-i+1};
                        end
                        P = lyap(A, B, C);
                        QLDPn{n+1} = P;
                        dqlMoms(n) = sum(inis*P);
                    end
                    dqlMoms = MomsFromFactorialMoms(dqlMoms);
                    % now calculate it at random time instance
                    pi = CTMCSolve (sD);
                    lambdak = sum(pi*D{k+1});
                    QLPn = {pi};
                    qlMoms = zeros(numOfQLMoms,1);
                    iTerm = inv(ones(N,1)*pi - sD);
                    for n=1:numOfQLMoms
                        sumP = sum(inis*QLDPn{n+1}) + n*(inis*QLDPn{n} - QLPn{n}*D{k+1}/lambdak)*iTerm*sum(D{k+1},2);
                        P = sumP*pi + n*(QLPn{n}*D{k+1} - inis*QLDPn{n}*lambdak)*iTerm;
                        QLPn{n+1} = P;
                        qlMoms(n) = sum(P);
                    end
                    qlMoms = MomsFromFactorialMoms(qlMoms);
                    res = qlMoms;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'qlDistr')
                    % DISTRIBUTION OF THE NUMBER OF JOBS
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLProbs = varargin{argIx+1};
                    sDk = D0;
                    for i=1:k-1
                        sDk = sDk + D{i+1};
                    end
                    % first calculate it at departure instants
                    Psid = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, sDk, 'P', precision);
                    Pn = {Psid};
                    dqlProbs = inis*Psid;
                    for n=1:numOfQLProbs-1
                        A = Qspp + Psid*Qsmp;
                        B = sDk + Qsmp*Psid;
                        C = Pn{n}*D{k+1};
                        for i=1:n-1
                            C = C + Pn{i+1}*Qsmp*Pn{n-i+1};
                        end
                        P = lyap(A, B, C);
                        Pn{n+1} = P;
                        dqlProbs = [dqlProbs; inis*P];
                    end
                    % now calculate it at random time instance
                    pi = CTMCSolve (sD);
                    lambdak = sum(pi*D{k+1});
                    iTerm = inv(-(sD-D{k+1}));
                    qlProbs = lambdak*dqlProbs(1,:)*iTerm;
                    for n=1:numOfQLProbs-1
                        P = (qlProbs(n,:)*D{k+1}+lambdak*(dqlProbs(n+1,:)-dqlProbs(n,:)))*iTerm;
                        qlProbs = [qlProbs; P];
                    end
                    qlProbs = sum(qlProbs,2);
                    res = qlProbs;
                    argIx = argIx + 1;
                else
                    error (['MMAPPH1PRPR: Unknown parameter ' varargin{argIx}])
                end
                if retIx>length(Ret)
                    Ret{retIx} = res;
                else
                    Ret{retIx} = [Ret{retIx}, res];
                end
                retIx = retIx+1;
                argIx = argIx + 1;
            end
        elseif k==K
            % step 3. calculate the performance measures
            % ==========================================   
            retIx = 1;
            argIx = 1;
            while argIx<=length(varargin)
                res = [];
                if any(ismember(eaten, argIx))
                    argIx = argIx + 1;
                    continue;
                elseif strcmp(varargin{argIx},'stMoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = varargin{argIx+1};
                    rtMoms = zeros(numOfSTMoms,1);
                    for i=1:numOfSTMoms
                        rtMoms(i) = factorial(i)*kappa*inv(-Kw)^(i+1)*sum(Bw,2);
                    end
                    res = rtMoms;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'stDistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = varargin{argIx+1};
                    rtDistr = [];
                    for t=stCdfPoints
                        rtDistr = [rtDistr; kappa*inv(-Kw)*(eye(size(Kw))-expm(Kw*t))*sum(Bw,2)];
                    end
                    res = rtDistr;
                    argIx = argIx + 1;
                elseif strcmp(varargin{argIx},'qlMoms') || strcmp(varargin{argIx},'qlDistr')
                    L = kron(sD-D{k+1},eye(M(k)))+kron(eye(N),S{k});
                    B = kron(eye(N),s{k}*sigma{k});
                    F = kron(D{k+1},eye(M(k)));
                    L0 = kron(sD-D{k+1},eye(M(k)));
                    R = QBDFundamentalMatrices (B, L, F, 'R', precision);
                    p0 = CTMCSolve(L0+R*B);
                    p0 = p0/sum(p0*inv(eye(size(R))-R));                    
                    if strcmp(varargin{argIx},'qlMoms')
                        % MOMENTS OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLMoms = varargin{argIx+1};
                        qlMoms = zeros(numOfQLMoms,1);
                        for i=1:numOfQLMoms
                            qlMoms(i) = sum(factorial(i)*p0*R^i*inv(eye(size(R))-R)^(i+1));
                        end
                        res = MomsFromFactorialMoms(qlMoms);
                    elseif strcmp(varargin{argIx},'qlDistr')        
                        % DISTRIBUTION OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLProbs = varargin{argIx+1};
                        qlProbs = p0;
                        for i=1:numOfQLProbs-1
                            qlProbs = [qlProbs; p0*R^i];
                        end
                        res = sum(qlProbs,2);                       
                    end
                    argIx = argIx + 1;
                else
                    error (['MMAPPH1PRPR: Unknown parameter ' varargin{argIx}])
                end                                        
                if retIx>length(Ret)
                    Ret{retIx} = res;
                else
                    Ret{retIx} = [Ret{retIx}, res];
                end
                retIx = retIx+1;
                argIx = argIx + 1;
            end
        end
    end
    varargout = Ret;
end