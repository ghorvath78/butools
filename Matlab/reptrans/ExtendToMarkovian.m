%  [beta, B] = ExtendToMarkovian(alpha, A, maxSize, precision)
%  
%  Assume we have an existing monocyclic (or acyclic) 
%  representation of a matrix-exponential distribution 
%  described by matrix A and vector alpha such that A is 
%  Markovian but alpha is not. This procedure appends an 
%  appropriate Erlang tail to the representation that makes 
%  the result Markovian (both the generator matrix and the 
%  initial vector parameter), while keeping the distribution 
%  the same. In [1]_ it is proven that this is always 
%  possible if the initial (alpha,A) defines a distribuion 
%  (non-negative density).
%  
%  Parameters
%  ----------
%  alpha : vector, shape (1,M)
%      The (non-Markovian) initial vector
%  A : matrix, shape (M,M)            
%      The (Markovian) transient generator.
%  maxSize : int, optional
%      The procedure stops if more than maxSize new phases 
%      are required. The default value is 100
%  precision : double, optional
%      The initial vector is considered to be valid if the 
%      smallest entry is greater than -precision. The
%      default value is 1e-14
%  
%  Returns
%  -------
%  beta : vector, shape (1,N)
%      The Markovian initial vector (N>=M)
%  B : matrix, shape (N,N)
%      The Markovian transient generator (N>=M).
%  
%  References
%  ----------
%  .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
%          of phase-type distributions," Stoch. Models 15, 
%          759-778 (1999)

function [beta, B] = ExtendToMarkovian (alpha, A, maxSize, precision)

    if ~exist('precision','var')
        precision = 1e-14;
    end

    if ~exist('maxSize','var')
        maxSize = 100;
    end

    function E = addErlangTail (D, len, mu)
        DN = size(D,1);
        E = zeros(DN+len);
        E(1:DN,1:DN) = D;
        E(DN,DN+1) = -sum(D(DN,:));
        for ei=1:len
            E(DN+ei,DN+ei) = -mu;
            if ei<len
                E(DN+ei,DN+ei+1) = mu;
            end
        end
    end

    function beta = inivecwithtail (gamma, G, tailLengthx, mux)
        vlen = size(G,1)+tailLengthx;
        beta = zeros(1, vlen);
        WG = eye(size(G))+G/mux;
        opv = gamma;
        clv = -sum(G/mux,2);
        for k=vlen:-1:size(G,1)+1
            beta(k) = opv*clv;
            opv = opv * WG;
        end
        beta(1:size(G,1)) = opv;
    end

    % initial value of t0upper
    t0lower = 0;
    t0upper = 1;
    beta = alpha * expm(A*t0upper);
    while min(beta)<-precision
        t0upper = t0upper * 2;
        beta = alpha * expm(A*t0upper);
    end

    % interval bisectioning
    while (t0upper - t0lower)/(t0upper + t0lower) > precision
        t0 = (t0upper + t0lower) / 2;
        beta = alpha * expm(A*t0);
        if min(beta)<-precision
            t0lower = t0;
        else
            t0upper = t0;
        end    
    end
    t0 = t0upper;

    % find optimal length and rate parameters of the Erlang tail

    % points towards t0 sometimes give fewer states
    % thus we try to increase t0 gradually till the number of states
    % decreases
    
    increment = 1.1;
    bestT0 = -1;
    bestLupper = -1;
    for i=1:100
        % initial value of Lupper
        Llower = 1;
        Lupper = 1;
        beta = inivecwithtail (alpha, A, Lupper, Lupper/t0);
        while min(beta)<-precision && Lupper<maxSize
            Lupper = Lupper * 2;
            beta = inivecwithtail (alpha, A, Lupper, Lupper/t0);
        end

        success =  min(beta)>=-precision;
        if success
            % interval bisectioning
            while Lupper - Llower > 1
                L = round((Lupper + Llower) / 2);
                beta = inivecwithtail (alpha, A, L, L/t0);
                if min(beta)<-precision
                    Llower = L;
                else
                    Lupper = L;
                end    
            end
        end
        
        if success 
            if bestLupper>=0 && Lupper>bestLupper % there was a successfull attempt before, and this one is worse
                break;                            % stop and keep what we have
            else                                  % otherwise go on and increment t0
                bestLupper = Lupper;
                bestT0 = t0;
                t0 = t0 * increment;                  
            end
        else
            if bestLupper>=0
                break;
            else
                t0 = t0 * increment;                  
            end
        end
    end
    t0 = bestT0;
    Lupper = bestLupper;
    
    if Lupper<0
        error ('No positive representation found up to the given size!');
    end
        
    % final result
    beta = inivecwithtail (alpha, A, Lupper, Lupper/t0);
    B = addErlangTail (A, Lupper, Lupper/t0);
end
