%  Ret = FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, ...)
%  
%  Returns various performane measures of a fluid queue
%  with independent fluid arrival and service processes.
%  
%  Two types of boundary behavior is available. If 
%  srv0stop=false, the output process evolves continuously
%  even if the queue is empty. If srv0stop=true, the 
%  output process slows down if there is fewer fluid in
%  the queue than it can serve. If the queue is empty
%  and the fluid input rate is zero, the output process
%  freezes till fluid arrives.
%  
%  Parameters
%  ----------
%  Qin : matrix, shape (N,N)
%      The generator of the background Markov chain 
%      corresponding to the input process
%  Rin : matrix, shape (N,N)
%      Diagonal matrix containing the fluid input rates
%      associated to the states of the input background 
%      process
%  Qout : matrix, shape (N,N)
%      The generator of the background Markov chain 
%      corresponding to the output process
%  Rout : matrix, shape (N,N)
%      Diagonal matrix containing the fluid output rates
%      associated to the states of the input background 
%      process
%  srv0stop : bool
%      If true, the service output process slows down if
%      there is fewer fluid in the queue than it can 
%      serve. If false, the output process evolves 
%      continuously.
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
%  
%  References
%  ----------
%  .. [1] Horvath G, Telek M, "Sojourn times in fluid queues 
%         with independent and dependent input and output 
%         processes PERFORMANCE EVALUATION 79: pp. 160-181, 2014.

function varargout = FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, varargin)

    % parse options
    prec = 1e-14;
    needST = 0;
    needQL = 0;
    eaten = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'prec')
            prec = varargin{i+1};
            eaten = [eaten, i, i+1];
        elseif length(varargin{i})>2 && strcmp(varargin{i}(1:2),'st')
            needST = 1;
        elseif length(varargin{i})>2 && strcmp(varargin{i}(1:2),'fl')
            needQL = 1;
        end
    end
    
    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   
    global BuToolsCheckPrecision;

    if BuToolsCheckInput && ~CheckGenerator(Qin,false)
        error('FluFlu: Generator matrix Qin is not Markovian!');
    end
    if BuToolsCheckInput && ~CheckGenerator(Qout,false)
        error('FluFlu: Generator matrix Qout is not Markovian!');
    end
    if BuToolsCheckInput && (any(diag(Rin)<-BuToolsCheckPrecision) || any(diag(Rout)<-BuToolsCheckPrecision))
        error('FluFlu: Fluid rates Rin and Rout must be non-negative !');
    end  
    
    Iin = eye(size(Qin));
    Iout = eye(size(Qout));

    if needQL
        Q = kron(Qin,Iout)+kron(Iin,Qout);
        if srv0stop
            Q0 = kron(Qin,Iout)+kron(Rin, pinv(Rout)*Qout);
        else
            Q0 = Q;
        end
        [mass0, ini, K, clo] = GeneralFluidSolve (Q, kron(Rin,Iout)-kron(Iin,Rout), Q0, prec);
    end
    if needST
        Rh = kron(Rin,Iout) - kron(Iin,Rout);
        Qh = kron(Qin, Rout) + kron(Rin, Qout);       
        [massh, inih, Kh, cloh] = GeneralFluidSolve (Qh, Rh, [], prec);

        % sojourn time density in case of 
        % srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda
        % srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu    

        lambda = sum(CTMCSolve(Qin)*Rin);
        mu = sum(CTMCSolve(Qout)*Rout);
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
            B = SimilarityMatrixForVectors(inv(-K)*sum(clo,2), ones(size(K,1),2));
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
            % convert result to PH representation
            Delta = diag(linsolve(Kh',-inih')); % Delta = diag (inih*inv(-Kh));
            A = inv(Delta)*Kh'*Delta;       
            if ~srv0stop        
                alpha = sum(Delta*cloh*kron(Rin,Iout)/lambda,2)';
            else
                alpha = sum(Delta*cloh*kron(Rin,Rout)/lambda/mu,2)';
            end        
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stDistrME')
            % convert result to ME representation
            if ~srv0stop
                B = SimilarityMatrixForVectors(sum(cloh*kron(Rin,Iout)/lambda,2), ones(size(Kh,1),1));
            else
                B = SimilarityMatrixForVectors(sum(cloh*kron(Rin,Rout)/lambda/mu,2), ones(size(Kh,1),1));
            end
            iB = inv(B);
            A = B*Kh*iB;
            alpha = inih*inv(-Kh)*iB;
            Ret{end+1} = alpha;
            Ret{end+1} = A;
        elseif strcmp(varargin{argIx},'stMoms')
            numOfMoms = varargin{argIx+1};
            argIx = argIx + 1;
            moms = zeros(1,numOfMoms);
            if srv0stop
                kclo = cloh*kron(Rin,Rout)/lambda/mu;
            else
                kclo = cloh*kron(Rin,Iout)/lambda;
            end
            iKh = inv(-Kh);
            for m=1:numOfMoms
                moms(m) = factorial(m)*sum(inih*iKh^(m+1)*kclo);
            end
            Ret{end+1} = moms;
        elseif strcmp(varargin{argIx},'stDistr')
            points = varargin{argIx+1};
            argIx = argIx + 1;
            values = zeros(size(points));
            if srv0stop
                kclo = cloh*kron(Rin,Rout)/lambda/mu;
            else
                kclo = cloh*kron(Rin,Iout)/lambda;
            end
            iKh = inv(-Kh);
            for p=1:length(points)
                values(p) = 1-sum(inih*expm(Kh*points(p))*iKh*kclo);
            end
            Ret{end+1} = values;
        else
            error (['FluFluQueue: Unknown parameter ' varargin{argIx}])
        end
        argIx = argIx + 1;
    end
    varargout = Ret;
end

