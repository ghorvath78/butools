%  Ret = FluidQueue(Q, Rin, Rout, ...)
%  
%  Returns various performane measures of a fluid queue.
%  
%  In a fluid queue there is a background continuous time
%  Markov chain (given by generator Q), and diagonal
%  matrix Rin (Rout) whose ith entry provides the 
%  fluid rate at which fluid enters the queue (can be 
%  served) while the background process is in state i.
%  
%  Parameters
%  ----------
%  Q : matrix, shape (N,N)
%      The generator of the background Markov chain
%  Rin : matrix, shape (N,N)
%      Diagonal matrix containing the fluid input rates
%      associated to the states of the background process
%  Rout : matrix, shape (N,N)
%      Diagonal matrix containing the fluid output rates
%      associated to the states of the background process
%  further parameters : 
%      The rest of the function parameters specify the options
%      and the performance measures to be computed.
%  
%      The supported performance measures and options in this 
%      function are:
%  
%      +----------------+--------------------+--------------------------------------+
%      | Parameter name | Input parameters   | Output                               |
%      +================+====================+======================================+
%      | "flMoms"       | Number of moments  | The moments of the fluid level       |
%      +----------------+--------------------+--------------------------------------+
%      | "flDistr"      | A vector of points | The fluid level distribution at      |
%      |                |                    | the requested points (cdf)           |
%      +----------------+--------------------+--------------------------------------+
%      | "flDistrME"    | None               | The vector-matrix parameters of the  |
%      |                |                    | matrix-exponentially distributed     |
%      |                |                    | fluid level distribution             |
%      +----------------+--------------------+--------------------------------------+
%      | "flDistrPH"    | None               | The vector-matrix parameters of the  |
%      |                |                    | matrix-exponentially distributed     |
%      |                |                    | fluid level distribution, converted  |
%      |                |                    | to a PH representation               |
%      +----------------+--------------------+--------------------------------------+
%      | "stMoms"       | Number of moments  | The sojourn time moments of fluid    |
%      |                |                    | drops                                |
%      +----------------+--------------------+--------------------------------------+
%      | "stDistr"      | A vector of points | The sojourn time distribution at the |
%      |                |                    | requested points (cummulative, cdf)  |
%      +----------------+--------------------+--------------------------------------+
%      | "stDistrME"    | None               | The vector-matrix parameters of the  |
%      |                |                    | matrix-exponentially distributed     |
%      |                |                    | sojourn time distribution            |
%      +----------------+--------------------+--------------------------------------+
%      | "stDistrPH"    | None               | The vector-matrix parameters of the  |
%      |                |                    | matrix-exponentially distributed     |
%      |                |                    | sojourn time distribution, converted |
%      |                |                    | to a PH representation               |
%      +----------------+--------------------+--------------------------------------+
%      | "prec"         | The precision      | Numerical precision to check if the  |
%      |                |                    | input is valid and it is also used   |
%      |                |                    | as a stopping condition when solving |
%      |                |                    | the Riccati equation                 |
%      +----------------+--------------------+--------------------------------------+
%      | "Q0"           | Matrix, shape(N,N) | The generator of the background      |
%      |                |                    | Markov chain when the fluid level is |
%      |                |                    | zero. If not given, Q0=Q is assumed  |
%      +----------------+--------------------+--------------------------------------+
%  
%  Returns
%  -------
%  Ret : list of the performance measures
%      Each entry of the list corresponds to a performance 
%      measure requested. If there is just a single item, 
%      then it is not put into a list.
%  
%  Notes
%  -----
%  "flDistrME" and "stDistrME" behave much better numerically than 
%  "flDistrPH" and "stDistrPH".

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

    global BuToolsCheckPrecision;
    
    if BuToolsCheckInput && ~CheckGenerator(Q,false)
        error('FluidQueue: Generator matrix Q is not Markovian!');
    end
    if BuToolsCheckInput && ~isempty(Q0) && ~CheckGenerator(Q0,false)
        error('FluidQueue: Generator matrix Q0 is not Markovian!');
    end
    if BuToolsCheckInput && (any(diag(Rin)<-BuToolsCheckPrecision) || any(diag(Rout)<-BuToolsCheckPrecision))
        error('FluidQueue: Fluid rates Rin and Rout must be non-negative !');
    end  

    [mass0, ini, K, clo] = GeneralFluidSolve (Q, Rin-Rout, Q0, prec);
    if needST
        N = size(Q,1);
        iniKi = linsolve(K',-ini')'; % iniki = ini*inv(-K);
        lambda = sum(mass0*Rin + iniKi*clo*Rin);
    end
    
    Ret = {};
    argIx = 1;
    while argIx<=length(varargin)
        if any(ismember(eaten, argIx))
            argIx = argIx + 1;
            continue;
        elseif strcmp(varargin{argIx},'flDistrPH')
            % transform it to PH
            Delta = diag(linsolve(K',-ini')); % Delta = diag (ini*inv(-K));
            A = inv(Delta)*K'*Delta;
            alpha = sum(clo,2)'*Delta;
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'flDistrME')
            % transform it to ME
            B = SimilarityMatrixForVectors(inv(-K)*sum(clo,2), ones(size(K,1),1));
            Bi = inv(B);
            alpha = ini*Bi;
            A = B*K*Bi;
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'flMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iK = inv(-K);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(ini*iK^(m+1)*clo);
            end
            Ret{end+1} = moms;
        elseif strcmp(varargin{argIx},'flDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            iK = inv(-K);
            for p=1:length(points)
                values(p) = sum(mass0) + sum(ini*(eye(size(K,1))-expm(K*points(p)))*iK*clo);
            end
            Ret{end+1} = values;
        elseif strcmp(varargin{argIx},'stDistrPH')
            % transform it to PH
            Delta = diag(iniKi/lambda);
            alpha = reshape(clo*Rin,1,N*length(ini))*kron(eye(N),Delta);
            A = kron(Rout, inv(Delta)*K'*Delta) + kron(Q, eye(size(K)));
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stDistrME')
            B = SimilarityMatrixForVectors(reshape(inv(-K)*clo*Rin,N*length(ini),1),ones(N*length(ini),1));
            Bi = inv(B);
            alpha = kron(ones(1,N), ini/lambda)*Bi;
            A = B*(kron(sparse(Q'),speye(size(K))) + kron(sparse(Rout),sparse(K)))*Bi;        
            Ret{end+1} = alpha;
            Ret{end+1} = A;
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
            Ret{end+1} = moms;
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
            Ret{end+1} = values;
        else
            error (['FluidQueue: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    varargout = Ret;
end

