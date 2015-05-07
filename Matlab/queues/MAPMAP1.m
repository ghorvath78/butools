function varargout = MAPMAP1(D0, D1, S0, S1, varargin)

    % parse options
    prec = 1e-14;
    needST = 0;
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'prec')
            prec = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif length(varargin{i})>2 && strcmp(varargin{i}(1:2),'st')
            needST = 1;
        end
    end
    
    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1,prec)
        error('MAPMAP1: The arrival process (D0,D1) is not a valid MAP representation!');
    end
    
    if BuToolsCheckInput && ~CheckMAPRepresentation(S0,S1,prec)
        error('MAPMAP1: The service process (S0,S1) is not a valid MAP representation!');
    end

    IA = eye(size(D0,1));
    IS = eye(size(S0,1));
    
    B = kron(IA,S1);
    L = kron(D0,IS)+kron(IA,S0);
    F = kron(D1,IS);
    L0 = kron(D0,IS);
    
    [pi0, R] = QBDSolve (B, L, F, L0, prec);
    N = length(pi0);
    I = eye(N);
    
    if needST
        % calculate the distribution of the age at departures
        U = L + R*B;
        Rh = inv(-U)*F;
        T = kron(IA,S0) + Rh * B;
        eta = pi0*F*inv(I-Rh);
        eta = eta/sum(eta);
    end
    
    Ret = {};
    retIx = 1;
    argIx = 1;
    while argIx<=length(varargin)
        if any(ismember(eaten, argIx))
            argIx = argIx + 1;
            continue;
        elseif strcmp(varargin{argIx},'qlDistrDPH')
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
            % transform it to PH representation
            beta = CTMCSolve(S0+S1, prec);
            theta = DTMCSolve(inv(-D0)*D1, prec);
            vv = kron(theta,beta);
            ix = 1:N;
            nz = ix(vv>prec);
            delta = diag(vv(nz));   
            alpha = ones(1,N)*B(nz,:)'*delta / sum(beta*S1);
            A = inv(delta)*T(nz,nz)'*delta;
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistrME')
            Ret{retIx} = {eta, T};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iT = inv(-T);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(eta*iT^m);
            end
            Ret{retIx} = moms;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            for p=1:length(points)
                values(p) = 1-sum(eta*expm(T*points(p)));
            end
            Ret{retIx} = values;
            retIx = retIx + 1;
        else
            error (['MAPMAP1: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    if length(Ret)==1 && iscell(Ret{1})
        varargout = Ret{1};
    else
        varargout = Ret;
    end
end
