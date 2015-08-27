%  Ret = MAPMAP1(D0, D1, S0, S1, ...)
%  
%  Returns various performane measures of a continuous time
%  MAP/MAP/1 queue.
%  
%  In a MAP/MAP/1 queue both the arrival and the service
%  processes are characterized by Markovian arrival 
%  processes.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape(N,N)
%      The transitions of the arrival MAP not accompanied by
%      job arrivals
%  D1 : matrix, shape(N,N)
%      The transitions of the arrival MAP accompanied by
%      job arrivals
%  S0 : matrix, shape(N,N)
%      The transitions of the service MAP not accompanied by
%      job service
%  S1 : matrix, shape(N,N)
%      The transitions of the service MAP accompanied by
%      job service
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
%  %  Notes
%  -----
%  "ncDistrMG" and "stDistrME" behave much better numerically than 
%  "ncDistrDPH" and "stDistrPH".

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

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1)
        error('MAPMAP1: The arrival process (D0,D1) is not a valid MAP representation!');
    end
    
    if BuToolsCheckInput && ~CheckMAPRepresentation(S0,S1)
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
            % transform it to PH representation
            beta = CTMCSolve(S0+S1);
            theta = DTMCSolve(inv(-D0)*D1);
            vv = kron(theta,beta);
            ix = 1:N;
            nz = ix(vv>prec);
            delta = diag(vv(nz));   
            alpha = ones(1,N)*B(nz,:)'*delta / sum(beta*S1);
            A = inv(delta)*T(nz,nz)'*delta;
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stDistrME')
            Ret{end+1} = eta;
            Ret{end+1} = T;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            iT = inv(-T);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(eta*iT^m);
            end
            Ret{end+1} = moms;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            for p=1:length(points)
                values(p) = 1-sum(eta*expm(T*points(p)));
            end
            Ret{end+1} = values;
        else
            error (['MAPMAP1: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    varargout = Ret;
end
