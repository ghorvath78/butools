function Ret = QBDQueue(B, L, F, L0, varargin)

    % parse options
    prec = 1e-15;
    needST = 0;
    for i=1:length(varargin)
        if strcmp(varargin{i},'prec')
            prec = varargin{i+1};
        elseif length(varargin{i})>2 && strcmp(varargin{i}(1:2),'st')
            needST = 1;
        end
    end
    
    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckGenerator(B+L+F,prec)
        error('QBDQueue: The matrix sum (B+L+F) is not a valid generator of a Markov chain!');
    end
    
    if BuToolsCheckInput && ~CheckGenerator(L0+F,prec)
        error('QBDQueue: The matrix sum (L0+F) is not a valid generator of a Markov chain!');
    end

    [pi0, R] = QBDSolve (B, L, F, L0, prec);
    N = length(pi0);
    I = eye(N);
    
    if needST
        U = L + R*B;
        Rh = inv(-U)*F;
        eta = pi0*F*inv(I-Rh);
        eta = eta/sum(eta);
        z = reshape(I,N*N,1);
    end
    
    Ret = {};
    retIx = 1;
    argIx = 1;
    while argIx<=length(varargin)
        if strcmp(varargin{argIx},'qlDistrDPH')
            % transform it to DPH
            alpha = pi0*R*inv(eye(N)-R);
            A = inv(diag(alpha))*R'*diag(alpha);           
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlDistrMG')
            % transform it to MG
            B = TransformToOnes(sum(inv(I-R)*R,2));
            Bi = inv(B);
            A = B*R*Bi;
            alpha = pi0*Bi;       
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iR = inv(I-R);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(pi0*iR^(m+1)*R^m);
            end
            Ret{retIx} = MomsFromFactorialMoms(moms);
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            for p=1:length(points)
                values(p) = sum(pi0*R^points(p));
            end
            Ret{retIx} = values;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistrPH')
            % transform to ph distribution
            ix = (1:N);
            nz = ix(eta>prec);
            Delta = diag(eta);
            A = kron(L+F,I(nz,nz)) + kron(B,inv(Delta(nz,nz))*Rh(nz,nz)'*Delta(nz,nz));
            alpha = z'*kron(I,Delta(:,nz));
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistrME')
            % transform it such that the closing vector is a vector of ones
            % this is the way butools accepts ME distributions
            Bm = TransformToOnes(z);
            Bmi = inv(Bm);
            A = Bm * (kron(L'+F',I) + kron(B',Rh)) * Bmi;
            alpha = kron(ones(1,N), eta) * Bmi;
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            Z = kron(L'+F',I)+kron(B',Rh);
            iZ = inv(-Z);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(kron(ones(1,N), eta)*iZ^(m+1)*(-Z)*z);
            end
            Ret{retIx} = moms;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            Z = kron(L'+F',I)+kron(B',Rh);
            for p=1:length(points)
                values(p) = 1-sum(kron(ones(1,N), eta)*expm(Z*points(p))*z);
            end
            Ret{retIx} = values;
            retIx = retIx + 1;
        else
            error (['QBDQueue: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    if length(Ret)==1
        Ret = Ret{1};
    end
end

