% out = MMAPPH1NPPR(D, sigma, S, perfmeas1, param1, perfmeas2, param2, ...)
%
% Solution of the MMAP[K]/PH[K]/1 non-preemptive priority queue
%
% Parameters
% ----------
% D : cell of matrices of shape (N,N), length (K+1)
%     The D0...DK matrices of the arrival process.
%     D1 corresponds to the lowest, DK to the highest priority.
% sigma : cell of row vectors, length (K)
%     The cell containing the initial probability vectors of the service
%     time distributions of the various customer types. The length of the
%     vectors does not have to be the same.
% S : cell of square matrices, length (K)
%     The transient generators of the phase type distributions representing
%     the service time of the jobs belonging to various types.
% list of perfmeas -- param pairs, where
%     perfmeas describes a performance measure to calculate, and the
%     corresponding param is the parameter of the calculation.
%     Supported performance measures are:
%     perfmeas = 'stmoms': returns the moments of the sojourn time. The
%     parameter is the number of moments to compute.
%     perfmeas = 'stdistr': returns the distribution of the sojourn time.
%     The parameter is the vector of points where the cdf is evaluated.
%     perfmeas = 'qlmoms': returns the moments of the number of jobs.
%     The parameter is the number of moments to compute.
%     perfmeas = 'qldistr': returns the distribution of the number of jobs.
%     The parameter is the upper bound.
% list of optionname -- value pairs, where the valid options are
%     optionname = 'erlMaxOrder': sets the maximal Erlang order used in the
%     erlangization procedure. The default value is 200.
%     optionname = 'precision': sets the numerical precision for the
%     solution of the Riccati and the matrix quadratic equations. The
%     default value is 1e-15.
%     optionname = 'classes': specifies which job classes are analyzed.
%     The default value is 1:K (all classes are analyzed).
%
% Returns
% -------
% Ret : cell of matrices
%     Each entry of the cell corresponds to a performance measure
%     requested.
%     Each entry is a matrix, where the columns belong to the various job
%     types.
%
function Ret = MMAPPH1NPPR(D, sigma, S, varargin)
    
    K = length(D)-1;

    % parse options
    erlMaxOrder = 200;
    precision = 1e-15;
    classes = 1:K;
    for i=1:length(varargin)
        if strcmp(varargin{i},'erlMaxOrder')
            erlMaxOrder = varargin{i+1};
        elseif strcmp(varargin{i},'precision')
            precision = varargin{i+1};
        elseif strcmp(varargin{i},'classes')
            classes = varargin{i+1};
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

    % step 1. solution of the workload process of the joint queue
    % ===========================================================
    Qwmm = D0;
    Qwpp = zeros(N*sum(M));
    Qwmp = zeros(N,N*sum(M));
    Qwpm = zeros(N*sum(M),N);
    kix = 1;
    for i=1:K
        bs = N*M(i);
        Qwpp(kix:kix+bs-1,kix:kix+bs-1) = kron(eye(N),S{i});
        Qwmp(:,kix:kix+bs-1) = kron(D{i+1},sigma{i});
        Qwpm(kix:kix+bs-1,:) = kron(eye(N),s{i});
        kix = kix + bs;
    end

    % calculate fundamental matrices
    Psiw = SolveRiccati (Qwpp, Qwpm, Qwmp, Qwmm, precision);
    Kw = Qwpp + Psiw*Qwmp;

    % calculate boundary vector
    Uw = Qwmm + Qwmp*Psiw;
    Ua = ones(N,1) + 2*sum(Qwmp*inv(-Kw),2);   
    pm = linsolve ([Uw,Ua]', [zeros(1,N),1]')';       
          
    ro =  ((1-sum(pm))/2)/(sum(pm)+(1-sum(pm))/2); % calc idle time with weight=1, and the busy time with weight=1/2
    kappa = pm/sum(pm);
    
    pi = CTMCSolve (sD);
    lambda = zeros(K,1);
    for i=1:K
        lambda(i) = sum(pi*D{i+1});
    end
    
    Psiw = {};
    Qwmp = {};  Qwzp = {};  Qwpp = {};
    Qwmz = {};  Qwpz = {};  Qwzz = {};
    Qwmm = {};  Qwpm = {};  Qwzm = {};
    for k=1:K
        % step 2. construct a workload process for classes k...K
        % ======================================================
        Mlo = sum(M(1:k-1));
        Mhi = sum(M(k:K));

        Qkwpp = zeros(N*Mlo*Mhi+N*Mhi);
        Qkwpz = zeros(N*Mlo*Mhi+N*Mhi, N*Mlo); 
        Qkwpm = zeros(N*Mlo*Mhi+N*Mhi, N);
        Qkwmz = zeros(N, N*Mlo);
        Qkwmp = zeros(N, N*Mlo*Mhi+N*Mhi);
        Dlo = D0;
        for i=1:k-1
            Dlo = Dlo + D{i+1};
        end
        Qkwmm = Dlo;
        Qkwzp = zeros(N*Mlo, N*Mlo*Mhi+N*Mhi);
        Qkwzm = zeros(N*Mlo, N);
        Qkwzz = zeros(N*Mlo);
        kix = 1;
        for i=k:K
            kix2 = 1;
            for j=1:k-1
                bs = N*M(j)*M(i);
                bs2 = N*M(j);
                Qkwpp(kix:kix+bs-1,kix:kix+bs-1) = kron(eye(N),kron(eye(M(j)),S{i}));
                Qkwpz(kix:kix+bs-1,kix2:kix2+bs2-1) = kron(eye(N),kron(eye(M(j)),s{i}));
                Qkwzp(kix2:kix2+bs2-1,kix:kix+bs-1) = kron(D{i+1},kron(eye(M(j)), sigma{i}));
                kix = kix + bs;
                kix2 = kix2 + bs2;
            end
        end
        for i=k:K
            bs = N*M(i);
            Qkwpp(kix:kix+bs-1,kix:kix+bs-1) = kron(eye(N),S{i});
            Qkwpm(kix:kix+bs-1,:) = kron(eye(N),s{i});
            Qkwmp(:,kix:kix+bs-1) = kron(D{i+1},sigma{i});
            kix = kix + bs;
        end
        kix = 1;
        for j=1:k-1
            bs = N*M(j);
            Qkwzz(kix:kix+bs-1,kix:kix+bs-1) = kron(Dlo, eye(M(j))) + kron(eye(N), S{j});
            Qkwzm(kix:kix+bs-1,:) = kron(eye(N), s{j});
            kix = kix + bs;
        end

        Psikw = SolveRiccati (Qkwpp+Qkwpz*inv(-Qkwzz)*Qkwzp, Qkwpm+Qkwpz*inv(-Qkwzz)*Qkwzm, Qkwmp, Qkwmm, precision);
        Psiw{k} = Psikw;
        
        Qwzp{k} = Qkwzp;    Qwmp{k} = Qkwmp;    Qwpp{k} = Qkwpp;
        Qwmz{k} = Qkwmz;    Qwpz{k} = Qkwpz;    Qwzz{k} = Qkwzz;
        Qwmm{k} = Qkwmm;    Qwpm{k} = Qkwpm;    Qwzm{k} = Qkwzm;
    end
    
    % step 3. calculate Phi vectors
    % =============================
    lambdaS = sum(lambda);
    phi{1} = (1-ro)*kappa*(-D0) / lambdaS;
    q0{1} = [];
    qL{1} = [];
    for k=1:K-1
        sDk = D{1};
        for j=1:k
            sDk = sDk + D{j+1};
        end
        % pk
        pk(k) = sum(lambda(1:k))/lambdaS - (1-ro)*kappa*sum(sDk,2)/lambdaS;
        % A^(k,1)
        Qwzpk = Qwzp{k+1};
        vix = 1;
        for ii=1:k
            bs = N*M(ii);
            V1 = Qwzpk(vix:vix+bs-1,:);
            Ak{ii} = kron(I,sigma{ii}) * inv(-kron(sDk,eye(M(ii)))-kron(I,S{ii})) * (kron(I,s{ii}) + V1*Psiw{k+1});
            vix = vix+ bs;
        end
        % B^k
        Qwmpk = Qwmp{k+1};
        Bk = Qwmpk * Psiw{k+1};
        ztag = phi{1}*(inv(-D0)*D{k+1}*Ak{k} - Ak{1} + inv(-D0)*Bk);
        for i=1:k-1
            ztag = ztag + phi{i+1}*(Ak{i}-Ak{i+1}) + phi{1}*inv(-D0)*D{i+1}*Ak{i};
        end
        Mx = eye(size(Ak{k}))-Ak{k};
        Mx(:,1) = ones(N,1);
        phi{k+1} = [pk(k), ztag(2:end)]*inv(Mx);
        
        q0{k+1} = phi{1}*inv(-D0);
        qL{k+1} = [];
        for ii=1:k
            qLii = (phi{ii+1} - phi{ii} + phi{1}*inv(-D0)*D{ii+1}) * kron(I,sigma{ii}) * inv(-kron(sDk,eye(M(ii)))-kron(I,S{ii}));
            qL{k+1} = [qL{k+1}, qLii];
        end
    end
    
    % step 4. calculate performance measures
    % ======================================
    Ret = {};
    for k=classes
        sD0k = D0;
        for i=1:k-1
            sD0k = sD0k + D{i+1};
        end
        if k<K
            % step 4.1 calculate distribution of the workload process right 
            % before the arrivals of class k jobs
            % ============================================================
            Kw = Qwpp{k}+Qwpz{k}*inv(-Qwzz{k})*Qwzp{k} + Psiw{k}*Qwmp{k};                   
            BM = []; CM = []; DM = [];
            for i=1:k-1
                BM = blkdiag(BM,kron(I,S{i}));
                CM = [CM; kron(I,s{i})];
                DM = blkdiag(DM,kron(D{k+1},eye(M(i))));
            end
            Kwu = [Kw, (Qwpz{k}+Psiw{k}*Qwmz{k})*inv(-Qwzz{k})*DM; zeros(size(BM,1),size(Kw,2)), BM];
            Bwu = [Psiw{k}*D{k+1}; CM];
            if k>1
                iniw = [q0{k}*Qwmp{k}+qL{k}*Qwzp{k}, qL{k}*DM];
                pwu = q0{k}*D{k+1};
            else
                iniw = pm*Qwmp{k};
                pwu = pm*D{k+1};
            end
            norm = sum(pwu) + sum(iniw*inv(-Kwu)*Bwu);
            pwu = pwu / norm;
            iniw = iniw / norm;

            % step 4.2 create the fluid model whose first passage time equals the
            % WAITING time of the low prioroity customers
            % ==================================================================
            KN = size(Kwu,1);
            Qspp = zeros(KN+N*sum(M(k+1:end)));
            Qspm = zeros(KN+N*sum(M(k+1:end)), N);
            Qsmp = zeros(N, KN+N*sum(M(k+1:end)));
            Qsmm = sD0k + D{k+1};
            kix = 1;
            for i=k+1:K
                bs = N*M(i);
                Qspp(KN+kix:KN+kix+bs-1,KN+kix:KN+kix+bs-1) = kron(I,S{i});
                Qspm(KN+kix:KN+kix+bs-1,:) = kron(I,s{i});
                Qsmp(:,KN+kix:KN+kix+bs-1) = kron(D{i+1},sigma{i});
                kix = kix + bs;
            end
            Qspp(1:KN,1:KN) = Kwu;
            Qspm(1:KN,:) = Bwu;           
            inis = [iniw, zeros(1,N*sum(M(k+1:end)))];

            % calculate fundamental matrix
            Psis = SolveRiccati (Qspp, Qspm, Qsmp, Qsmm, precision);
            
            % step 4.3. calculate the performance measures
            % ==========================================   
            retIx = 1;
            for va=1:length(varargin)/2
                perfmeas = varargin{va*2-1};
                param = varargin{va*2};
                res = [];
                if strcmp(perfmeas,'stmoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = param;
                    % calculate waiting time moments
                    Pn = {Psis};
                    wtMoms = zeros(numOfSTMoms,1);
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
                        wtMoms(n) = sum(inis*P*(-1)^n) / 2^n;
                    end
                    % calculate RESPONSE time moments
                    Pnr = {sum(inis*Pn{1})*sigma{k}}; %{sigmaL}; 
                    rtMoms = zeros(numOfSTMoms,1);
                    for n=1:numOfSTMoms
                        P =  n*Pnr{n}*inv(-S{k}) + (-1)^n*sum(inis*Pn{n+1})*sigma{k} / 2^n;
                        Pnr{n+1} = P;
                        rtMoms(n) = sum(P)+sum(pwu)*factorial(n)*sum(sigma{k}*inv(-S{k})^n);
                    end
                    res = rtMoms;
                elseif strcmp(perfmeas,'stdistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = param;
                    for t=stCdfPoints
                        L = erlMaxOrder;
                        lambdae = L/t/2;
                        Psie = SolveRiccati (Qspp-lambdae*eye(size(Qspp)), Qspm, Qsmp, Qsmm-lambdae*eye(size(Qsmm)), precision);
                        Pn = {Psie};
                        pr = (sum(pwu) + sum(inis*Psie)) * (1-sum(sigma{k}*inv(eye(size(S{k}))-S{k}/2/lambdae)^L));
                        for n=1:L-1
                            A = Qspp + Psie*Qsmp - lambdae*eye(size(Qspp));
                            B = Qsmm + Qsmp*Psie - lambdae*eye(size(Qsmm));
                            C = 2*lambdae*Pn{n};
                            for i=1:n-1
                                C = C + Pn{i+1}*Qsmp*Pn{n-i+1};
                            end
                            P = lyap(A, B, C);
                            Pn{n+1} = P;
                            pr = pr + sum(inis*P) * (1-sum(sigma{k}*inv(eye(size(S{k}))-S{k}/2/lambdae)^(L-n)));
                        end
                        res = [res; pr];
                    end
                elseif strcmp(perfmeas,'qlmoms') || strcmp(perfmeas,'qldistr')
                    W = inv(-kron(sD-D{k+1},eye(M(k)))-kron(I,S{k}))*kron(D{k+1},eye(M(k)));
                    iW = inv(eye(size(W))-W);
                    w = kron(eye(N),sigma{k});
                    omega = inv(-kron(sD-D{k+1},eye(M(k)))-kron(I,S{k}))*kron(I,s{k});
                    if strcmp(perfmeas,'qlmoms')
                        % MOMENTS OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLMoms = param;
                        % calculate queue length moments at departures
                        Psii = {Psis};
                        QLDPn = {inis*Psii{1}*w*iW};
                        for n=1:numOfQLMoms
                            A = Qspp + Psis*Qsmp;
                            B = Qsmm + Qsmp*Psis;
                            C = n*Psii{n}*D{k+1};
                            bino = 1;
                            for i=1:n-1
                                bino = bino * (n-i+1) / i;
                                C = C + bino * Psii{i+1}*Qsmp*Psii{n-i+1};
                            end
                            P = lyap(A, B, C);
                            Psii{n+1} = P;
                            QLDPn{n+1} = n*QLDPn{n}*iW*W + inis*P*w*iW;
                        end
                        for n=0:numOfQLMoms
                            QLDPn{n+1} = (QLDPn{n+1} + pwu*w*iW^(n+1)*W^n)*omega;
                        end                    
                        % now calculate it at random time instance
                        QLPn = {pi};
                        qlMoms = zeros(numOfQLMoms,1);
                        iTerm = inv(ones(N,1)*pi - sD);
                        for n=1:numOfQLMoms
                            sumP = sum(QLDPn{n+1}) + n*(QLDPn{n} - QLPn{n}*D{k+1}/lambda(k))*iTerm*sum(D{k+1},2);
                            P = sumP*pi + n*(QLPn{n}*D{k+1} - QLDPn{n}*lambda(k))*iTerm;
                            QLPn{n+1} = P;
                            qlMoms(n) = sum(P);
                        end
                        qlMoms = MomsFromFactorialMoms(qlMoms);
                        res = qlMoms;
                    elseif strcmp(perfmeas,'qldistr')
                        % DISTRIBUTION OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLProbs = param;
                        Psid = SolveRiccati (Qspp, Qspm, Qsmp, sD0k);
                        Pn = {Psid};
                        XDn = inis*Psid*w;
                        dqlProbs = (XDn+pwu*w)*omega;
                        for n=1:numOfQLProbs-1
                            A = Qspp + Psid*Qsmp;
                            B = sD0k + Qsmp*Psid;
                            C = Pn{n}*D{k+1};
                            for i=1:n-1
                                C = C + Pn{i+1}*Qsmp*Pn{n-i+1};
                            end
                            P = lyap(A, B, C);
                            Pn{n+1} = P;
                            XDn = XDn*W + inis*P*w;
                            dqlProbs = [dqlProbs; (XDn+pwu*w*W^n)*omega];
                        end
                        iTerm = inv(-(sD-D{k+1}));
                        qlProbs = lambda(k)*dqlProbs(1,:)*iTerm;
                        for n=1:numOfQLProbs-1
                            P = (qlProbs(n,:)*D{k+1}+lambda(k)*(dqlProbs(n+1,:)-dqlProbs(n,:)))*iTerm;
                            qlProbs = [qlProbs; P];
                        end
                        qlProbs = sum(qlProbs,2);
                        res = qlProbs;
                    end
                end
                if retIx>length(Ret)
                    Ret{retIx} = res;
                else
                    Ret{retIx} = [Ret{retIx}, res];
                end
                retIx = retIx+1;
            end           
        elseif k==K
            % step 3. calculate the performance measures
            % ==========================================   
            retIx = 1;
            for va=1:length(varargin)/2
                perfmeas = varargin{va*2-1};
                param = varargin{va*2};
                res = [];
                if strcmp(perfmeas,'stmoms') || strcmp(perfmeas,'stdistr')
                    Kw = Qwpp{k}+Qwpz{k}*inv(-Qwzz{k})*Qwzp{k} + Psiw{k}*Qwmp{k};                   
                    AM = []; BM = []; CM = []; DM = [];
                    for i=1:k-1
                        AM = blkdiag(AM, kron(ones(N,1),kron(eye(M(i)),s{k})));
                        BM = blkdiag(BM,S{i});
                        CM = [CM; s{i}];
                        DM = blkdiag(DM,kron(D{k+1},eye(M(i))));
                    end
                    Z = [Kw, [AM;zeros(N*M(k),size(AM,2))]; zeros(size(BM,1),size(Kw,2)), BM];
                    z = [zeros(size(AM,1),1); kron(ones(N,1),s{k}); CM];
                    iniw = [q0{k}*Qwmp{k}+qL{k}*Qwzp{k}, zeros(1,size(BM,1))];
                    zeta = iniw/sum(iniw*inv(-Z)*z);   
                    if strcmp(perfmeas,'stmoms')
                        % MOMENTS OF THE SOJOURN TIME
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfSTMoms = param;
                        rtMomsH = zeros(numOfSTMoms,1);
                        for i=1:numOfSTMoms
                            rtMomsH(i) = factorial(i)*zeta*inv(-Z)^(i+1)*z;
                        end
                        res = rtMomsH;
                    elseif strcmp(perfmeas,'stdistr') 
                        % DISTRIBUTION OF THE SOJOURN TIME
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        stCdfPoints = param;
                        rtDistr = [];
                        for t=stCdfPoints
                            rtDistr = [rtDistr; zeta*inv(-Z)*(eye(size(Z))-expm(Z*t))*z];
                        end
                        res = rtDistr;
                    end
                elseif strcmp(perfmeas,'qlmoms') || strcmp(perfmeas,'qldistr')
                    L = zeros(N*sum(M));
                    B = zeros(N*sum(M));
                    F = zeros(N*sum(M));
                    kix = 1;
                    for i=1:K
                        bs = N*M(i);
                        F(kix:kix+bs-1,kix:kix+bs-1) = kron(D{k+1},eye(M(i)));
                        L(kix:kix+bs-1,kix:kix+bs-1) = kron(sD0k,eye(M(i))) + kron(I,S{i});
                        if i<K
                            L(kix:kix+bs-1,N*sum(M(1:k-1))+1:end) = kron(I,s{i}*sigma{k});
                        else
                            B(kix:kix+bs-1,N*sum(M(1:k-1))+1:end) = kron(I,s{i}*sigma{k});
                        end
                        kix = kix + bs;
                    end
                    R = SolveMatrixQuadratic (F, L, B, precision);
                    p0 = [qL{k}, q0{k}*kron(I,sigma{k})];
                    p0 = p0/sum(p0*inv(eye(size(R))-R));
                    if strcmp(perfmeas,'qlmoms')
                        % MOMENTS OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLMoms = param;
                        qlMoms = zeros(numOfQLMoms,1);
                        for i=1:numOfQLMoms
                            qlMoms(i) = sum(factorial(i)*p0*R^i*inv(eye(size(R))-R)^(i+1));
                        end
                        res = MomsFromFactorialMoms(qlMoms);
                    elseif strcmp(perfmeas,'qldistr')        
                        % DISTRIBUTION OF THE NUMBER OF JOBS
                        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLProbs = param;
                        qlProbs = p0;
                        for i=1:numOfQLProbs-1
                            qlProbs = [qlProbs; p0*R^i];
                        end
                        res = sum(qlProbs,2);                       
                    end
                end                                                       
                if retIx>length(Ret)
                    Ret{retIx} = res;
                else
                    Ret{retIx} = [Ret{retIx}, res];
                end
                retIx = retIx+1;
            end           
        end
    end   
end