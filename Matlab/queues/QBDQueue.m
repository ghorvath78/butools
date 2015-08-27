%  Ret = QBDQueue(B, L, F, L0, ...)
%  
%  Returns various performane measures of a continuous time
%  QBD queue.
%  
%  QBD queues have a background continuous time Markov chain
%  with generator Q whose the transitions can be partitioned
%  into three sets: transitions accompanied by an arrival
%  of a new job (F, forward), transitions accompanied by 
%  the service of the current job in the server (B, 
%  backward) and internal transitions (L, local). 
%  Thus we have Q=B+L+F. L0 is the matrix of local 
%  transition rates if the queue is empty.
%  
%  Parameters
%  ----------
%  B : matrix, shape(N,N)
%      Transitions of the background process accompanied by 
%      the service of the current job in the server
%  L : matrix, shape(N,N)
%      Internal transitions of the background process 
%      that do not generate neither arrival nor service
%  F : matrix, shape(N,N)
%      Transitions of the background process accompanied by 
%      an arrival of a new job
%  L0 : matrix, shape(N,N)
%      Internal transitions of the background process when
%      there are no jobs in the queue
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
%      | "ncDistrMG"    | None               | The vector-matrix parameters of the    |
%      |                |                    | matrix-geometric distribution of the   |
%      |                |                    | number of customers in the system      |
%      +----------------+--------------------+----------------------------------------+
%      | "ncDistrDPH"   | None               | The vector-matrix parameters of the    |
%      |                |                    | matrix-geometric distribution of the   |
%      |                |                    | number of customers in the system,     |
%      |                |                    | converted to a discrete PH             |
%      |                |                    | representation                         |
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
%      |                |                    | condition when solving the             |
%      |                |                    | matrix-quadratic equation              |
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
%      measure requested. If there is just a single item, 
%      then it is not put into a list.
%  
%  Notes
%  -----
%  "ncDistrMG" and "stDistrMG" behave much better numerically than 
%  "ncDistrDPH" and "stDistrPH".

function varargout = QBDQueue(B, L, F, L0, varargin)

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

    if BuToolsCheckInput && ~CheckGenerator(B+L+F)
        error('QBDQueue: The matrix sum (B+L+F) is not a valid generator of a Markov chain!');
    end
    
    if BuToolsCheckInput && ~CheckGenerator(L0+F)
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
    argIx = 1;
    while argIx<=length(varargin)
        if any(ismember(eaten, argIx))
            argIx = argIx + 1;
            continue;
        elseif strcmp(varargin{argIx},'ncDistrDPH')
            % transform it to DPH
            alpha = pi0*R*inv(eye(N)-R);
            A = inv(diag(alpha))*R'*diag(alpha);           
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'ncDistrMG')
            % transform it to MG
            B = SimilarityMatrixForVectors(sum(inv(I-R)*R,2), ones(N,1));
            Bi = inv(B);
            A = B*R*Bi;
            alpha = pi0*Bi;       
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'ncMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iR = inv(I-R);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(pi0*iR^(m+1)*R^m);
            end
            Ret{end+1} = MomsFromFactorialMoms(moms);
        elseif strcmp(varargin{argIx},'ncDistr')
            numOfQLProbs = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(1,numOfQLProbs);
            values(1) = sum(pi0);
            RPow = I;
            for p=1:numOfQLProbs-1
                RPow = RPow * R;
                values(p+1) = sum(pi0*RPow);
            end
            Ret{end+1} = values;
        elseif strcmp(varargin{argIx},'stDistrPH')
            % transform to ph distribution
            ix = (1:N);
            nz = ix(eta>prec);
            Delta = diag(eta);
            A = kron(L+F,I(nz,nz)) + kron(B,inv(Delta(nz,nz))*Rh(nz,nz)'*Delta(nz,nz));
            alpha = z'*kron(I,Delta(:,nz));
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stDistrME')
            % transform it such that the closing vector is a vector of ones
            % this is the way butools accepts ME distributions
            Bm = SimilarityMatrixForVectors(z,ones(length(z),1));
            Bmi = inv(Bm);
            A = Bm * (kron(L'+F',I) + kron(B',Rh)) * Bmi;
            alpha = kron(ones(1,N), eta) * Bmi;
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            Z = kron(L'+F',I)+kron(B',Rh);
            iZ = inv(-Z);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(kron(ones(1,N), eta)*iZ^(m+1)*(-Z)*z);
            end
            Ret{end+1} = moms;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            Z = kron(L'+F',I)+kron(B',Rh);
            for p=1:length(points)
                values(p) = 1-sum(kron(ones(1,N), eta)*expm(Z*points(p))*z);
            end
            Ret{end+1} = values;
        else
            error (['QBDQueue: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    varargout = Ret;
end

