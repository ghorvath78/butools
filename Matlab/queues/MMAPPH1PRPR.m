% out = MMAPPH1PRPR(D, sigma, S, perfmeas1, param1, perfmeas2, param2, ...)
%
% Solution of the MMAP[K]/PH[K]/1 preemptive resume priority queue
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
function Ret = MMAPPH1PRPR(D, sigma, S, varargin)
    
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
        Psiw = SolveRiccati (Qwpp, Qwpm, Qwmp, Qwmm, precision);
        Kw = Qwpp + Psiw*Qwmp;

        % calculate boundary vector
        Uw = Qwmm + Qwmp*Psiw;
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
            Psis = SolveRiccati (Qspp, Qspm, Qsmp, Qsmm, precision);            
            
            % step 3. calculate the performance measures
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
                elseif strcmp(perfmeas,'stdistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = param;
                    for t=stCdfPoints
                        L = erlMaxOrder;
                        lambda = L/t/2;
                        Psie = SolveRiccati (Qspp-lambda*eye(size(Qspp)), Qspm, Qsmp, Qsmm-lambda*eye(size(Qsmm)), precision);
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
                elseif strcmp(perfmeas,'qlmoms')
                    % MOMENTS OF THE NUMBER OF JOBS
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLMoms = param;
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
                elseif strcmp(perfmeas,'qldistr')
                    % DISTRIBUTION OF THE NUMBER OF JOBS
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLProbs = param;
                    sDk = D0;
                    for i=1:k-1
                        sDk = sDk + D{i+1};
                    end
                    % first calculate it at departure instants
                    Psid = SolveRiccati (Qspp, Qspm, Qsmp, sDk, precision);
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
                if strcmp(perfmeas,'stmoms') 
                    % MOMENTS OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = param;
                    rtMoms = zeros(numOfSTMoms,1);
                    for i=1:numOfSTMoms
                        rtMoms(i) = factorial(i)*kappa*inv(-Kw)^(i+1)*sum(Bw,2);
                    end
                    res = rtMoms;
                elseif strcmp(perfmeas,'stdistr') 
                    % DISTRIBUTION OF THE SOJOURN TIME
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = param;
                    rtDistr = [];
                    for t=stCdfPoints
                        rtDistr = [rtDistr; kappa*inv(-Kw)*(eye(size(Kw))-expm(Kw*t))*sum(Bw,2)];
                    end
                    res = rtDistr;
                elseif strcmp(perfmeas,'qlmoms') || strcmp(perfmeas,'qldistr')
                    L = kron(sD-D{k+1},eye(M(k)))+kron(eye(N),S{k});
                    B = kron(eye(N),s{k}*sigma{k});
                    F = kron(D{k+1},eye(M(k)));
                    L0 = kron(sD-D{k+1},eye(M(k)));
                    R = SolveMatrixQuadratic (F, L, B, precision);
                    p0 = CTMCSolve(L0+R*B);
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