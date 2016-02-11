%  M = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method)
%  
%  Returns the fundamental matrices corresponding to the
%  given canonical Markov fluid model. Matrices Psi, K and
%  U are returned depending on the "matrices" parameter.
%  The canonical Markov fluid model is defined by the 
%  matrix blocks of the generator of the background Markov
%  chain partitioned according to the sign of the 
%  associated fluid rates (i.e., there are "+" and "-" states).
%  
%  Parameters
%  ----------
%  Fpp : matrix, shape (Np,Np)
%      The matrix of transition rates between states 
%      having positive fluid rates
%  Fpm : matrix, shape (Np,Nm)
%      The matrix of transition rates where the source
%      state has a positive, the destination has a 
%      negative fluid rate associated.
%  Fpm : matrix, shape (Nm,Np)
%      The matrix of transition rates where the source
%      state has a negative, the destination has a 
%      positive fluid rate associated.
%  Fpp : matrix, shape (Nm,Nm)
%      The matrix of transition rates between states 
%      having negative fluid rates
%  matrices : string
%      Specifies which matrices are required. 'P' means 
%      that only matrix Psi is needed. 'UK' means that
%      matrices U and K are needed. Any combinations of
%      'P', 'K' and 'U' are allowed, in any order.
%  precision : double, optional
%      The matrices are computed iteratively up to this
%      precision. The default value is 1e-14
%  maxNumIt : int, optional
%      The maximal number of iterations. The default value
%      is 50.
%  method : {"CR", "ADDA", "SDA"}, optional
%      The method used to solve the algebraic Riccati
%      equation (CR: cyclic reduction, ADDA: alternating-
%      directional doubling algorithm, SDA: structured
%      doubling algorithm). The default is "CR".
%  
%  Returns
%  -------
%  M : list of matrices
%      The list of calculated matrices in the order as
%      requested in the 'matrices' parameter.
%  
%  Notes
%  -----
%  Thanks to Benny Van Houdt for the implementation of the
%  Riccati solvers.

function varargout = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method)

    if ~exist('precision','var')
        precision = 1e-14;
    end
    if ~exist('maxNumIt','var')
        maxNumIt = 150;
    end
    if ~exist('method','var')
        method = 'ADDA';
    end

    if size(Fpp,1)==0
        numit = 0;
        Psi = zeros(0,size(Fmm,1));
    elseif strcmp(method,'CR')
        F = [Fpp,Fpm;Fmp,Fmm];
        sp = size(Fpp,1);
        sm = size(Fmm,1);
        I = eye(sm+sp);
        Pf = I + F / max(-diag(F));
        A1 = [Pf(1:sp,1:sp)/2.0, zeros(sp,sm); Pf(sp+1:end,1:sp), zeros(sm,sm)];
        % Optimized cyclic reduction
        A1hat=A1;
        R1 = eye(sp) / 2.0;
        N = [Pf(1:sp,sp+1:end)/2.0; Pf(sp+1:sp+sm,sp+1:sp+sm)];
        numit = 0;
        while min(norm(R1,inf),norm(N,inf))> precision && numit<maxNumIt
            Xmat = inv(I-A1);
            B = Xmat*N;
            S1 = Xmat(sp+1:end,1:sp)*R1;
            U = R1*B(1:sp,:);
            A1 = A1 + [N*S1, [U; zeros(sm,sm)]];
            A1hat(1:sp,sp+1:end) = A1hat(1:sp,sp+1:end) + U;
            N = N * B(sp+1:end,:);
            R1 = R1 * Xmat(1:sp,1:sp) * R1;
            numit = numit + 1;
        end
        Ghat = inv(I-A1hat) * [Pf(1:sp,sp+1:end)/2.0; Pf(sp+1:end,sp+1:end)];
        Psi = Ghat(1:sp,:);
    elseif strcmp(method,'ADDA') || strcmp(method,'SDA')
        % via ADDA algorithm (Wang, Wang, Li 2011)
        A = -Fpp;
        B = Fpm;
        C = Fmp;
        D = -Fmm;
        gamma1 = max(diag(A));
        gamma2 = max(diag(D));
        if strcmp(method,'SDA')
            gamma1 = max(gamma1,gamma2);
            gamma2 = gamma1;
        end
        sA = size(A,1);
        sD = size(D,1);
        IA = eye(sA);
        ID = eye(sD);
        A = A + gamma2*IA;
        D = D + gamma1*ID;
        Dginv = inv(D);
        Vginv = inv(D-C*inv(A)*B);
        Wginv = inv(A-B*Dginv*C);
        Eg = ID - (gamma1+gamma2)*Vginv;
        Fg = IA - (gamma1+gamma2)*Wginv;
        Gg = (gamma1+gamma2)*Dginv*C*Wginv;
        Hg = (gamma1+gamma2)*Wginv*B*Dginv;
        diff = 1.0;
        numit = 0;
        while diff > precision && numit < maxNumIt
            Vginv = Eg*inv(ID-Gg*Hg);
            Wginv = Fg*inv(IA-Hg*Gg);
            Gg = Gg + Vginv*Gg*Fg;
            Hg = Hg + Wginv*Hg*Eg;
            Eg = Vginv*Eg;
            Fg = Wginv*Fg;
            neg = norm(Eg,1);
            nfg = norm(Fg,1);
            if strcmp(method,'ADDA')
                eta = sqrt(nfg/neg);
                Eg = Eg*eta;
                Fg = Fg/eta;
                diff = neg*nfg;
            else
                diff = min(neg,nfg);
            end
            numit = numit + 1;
        end
        Psi = Hg;
    end
    
    global BuToolsVerbose;
    if numit == maxNumIt && BuToolsVerbose==true
        fprintf('Maximum Number of Iterations reached');
    end
    
    if BuToolsVerbose==true
        res_norm = norm (Fpm+Fpp*Psi+Psi*Fmm+Psi*Fmp*Psi, inf);
        fprintf('Final Residual Error for Psi: ');
        disp(res_norm);
    end

    ret = cell(1,length(matrices));
    for i=1:length(matrices)
        if matrices(i)=='P'
            ret{i} = Psi;
        elseif matrices(i)=='K'
            ret{i} = Fpp+Psi*Fmp;
        elseif matrices(i)=='U'
            ret{i} = Fmm+Fmp*Psi;
        end
    end
    varargout=ret;
end

