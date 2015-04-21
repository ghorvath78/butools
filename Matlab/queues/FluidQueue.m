function varargout = FluidQueue(Q, Rin, Rout, varargin)

    % parse options
    prec = 1e-14;
    needST = 0;
    Q0 = [];
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'prec')
            prec = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif strcmp(varargin{i},'Q0')
            Q0 = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif length(varargin{i})>2 && strcmp(varargin{i}(1:2),'st')
            needST = 1;
        end
    end
    
    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckGenerator(Q,false,prec)
        error('FluidQueue: Generator matrix Q is not Markovian!');
    end
    if BuToolsCheckInput && ~isempty(Q0) && ~CheckGenerator(Q0,false,prec)
        error('FluidQueue: Generator matrix Q0 is not Markovian!');
    end
    if BuToolsCheckInput && (any(diag(Rin)<-prec) || any(diag(Rout)<-prec))
        error('FluidQueue: Fluid rates Rin and Rout must be non-negative !');
    end  

    [mass0, ini, K, clo] = GeneralFluidSolve (Q, Rin-Rout, Q0, prec);
    if needST
        N = size(Q,1);
        iniKi = linsolve(K',-ini')'; % iniki = ini*inv(-K);
        lambda = sum(mass0*Rin + iniKi*clo*Rin);
    end
    
    Ret = {};
    retIx = 1;
    argIx = 1;
    while argIx<=length(varargin)
        if any(ismember(eaten, argIx))
            argIx = argIx + 1;
            continue;
        elseif strcmp(varargin{argIx},'qlDistrPH')
            % transform it to PH
            Delta = diag(linsolve(K',-ini')); % Delta = diag (ini*inv(-K));
            A = inv(Delta)*K'*Delta;
            alpha = sum(clo,2)'*Delta;
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlDistrME')
            % transform it to ME
            B = TransformToOnes(inv(-K)*sum(clo,2));
            Bi = inv(B);
            alpha = ini*Bi;
            A = B*K*Bi;
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iK = inv(-K);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(ini*iK^(m+1)*clo);
            end
            Ret{retIx} = moms;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'qlDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            iK = inv(-K);
            for p=1:length(points)
                values(p) = sum(mass0) + sum(ini*(eye(size(K,1))-expm(K*points(p)))*iK*clo);
            end
            Ret{retIx} = values;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistrPH')
            % transform it to PH
            Delta = diag(iniKi/lambda);
            alpha = reshape(clo*Rin,1,N*length(ini))*kron(eye(N),Delta);
            A = kron(Rout, inv(Delta)*K'*Delta) + kron(Q, eye(size(K)));
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistrME')
            B = TransformToOnes(reshape(inv(-K)*clo*Rin,N*length(ini),1));
            Bi = inv(B);
            alpha = kron(ones(1,N), ini/lambda)*Bi;
            A = B*(kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K)))*Bi;        
            Ret{retIx} = {alpha, A};
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            Z = kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K));
            iZ = inv(-Z);
            kini = kron(ones(1,N), ini/lambda);
            kclo = reshape(inv(-K)*clo*Rin,N*length(ini),1);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(kini*iZ^(m+1)*(-Z)*kclo);
            end
            Ret{retIx} = moms;
            retIx = retIx + 1;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            Z = kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K));
            kini = kron(ones(1,N), ini/lambda);
            kclo = reshape(inv(-K)*clo*Rin,N*length(ini),1);
            for p=1:length(points)
                values(p) = 1-sum(kini*expm(Z*points(p))*kclo);
            end
            Ret{retIx} = values;
            retIx = retIx + 1;
        else
            error (['FluidQueue: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    if length(Ret)==1 && iscell(Ret{1})
        varargout = Ret{1};
    else
        varargout = Ret;
    end
end

