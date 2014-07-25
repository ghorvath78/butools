%  B = TransformToMonocyclic(A, maxSize, precision)
%  
%  Transforms an arbitrary matrix to a Markovian monocyclic 
%  matrix (see [1]_).
%  
%  Parameters
%  ----------
%  A : matrix, shape (N,N)
%      Matrix parameter of the initial representation
%  maxSize : int, optional
%      The maximal order of the resulting Markovian 
%      representation. The default value is 100
%  precision : double, optional
%      Matrix entries smaller than the precision are 
%      considered to be zeros. The default value is 1e-14
%      
%  Returns
%  -------
%  B : matrix, shape (M,M)
%      Transient generator matrix of the Markovian monocyclic
%      representation. Note that M>N if there are complex 
%      eigenvalues.
%  
%  Notes
%  -----    
%  Raises an exception if no Markovian monocyclic generator 
%  has been found up to the given size.
%  
%  References
%  ----------
%  .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
%          of phase-type distributions," Stoch. Models 15, 
%          759-778 (1999)

function B = TransformToMonocyclic (A, maxSize, precision)

    if ~exist('precision','var')
        precision = 1e-14;
    end

    if ~exist('maxSize','var')
        maxSize = 100;
    end

    function [evalues, repeats] = eigvalc(A)
        % eigval  Eigenvalues and their algebraic multiplicity.
        %
        % evalues = eigval(A) returns the distinct eigenvalues of A,
        % with duplicates removed and sorted in increasing order.
        %
        % [evalues, repeats] = eigval(A) also returns the row vector
        % repeats that gives the multiplicity of each eigenvalue.
        % The sum of the multiplicities is n.
        %
        % Examples: Let A = eye(n) and B = diag([3 4]).
        % For A, evalues is 1 and repeats is n.
        % For B, evalues is [4; 3]  and repeats is [1 1].

        tol = sqrt(eps);
        lambda = sort(eig(A));
        lambda = round(lambda/tol) * tol;
        %
        % lambda gives all n eigenvalues (repetitions included).
        %
        evalues = unique(lambda);
        n = length(lambda);
        d = length(evalues);
        A = ones(n, 1) * evalues';
        B = lambda * ones(1, d);
        MATCH = abs(A-B) <= tol;
        %
        % MATCH is an n by d zero matrix except
        % MATCH(i,j) = 1 when lambda(i) = evalues(j).
        % Summing the columns gives the row vector repeats.
        %
        repeats = sum(MATCH);
    end

    function febs = generateFEBs (A, evalues, repeats, maxSize, precision)
        % calculate the parameters of the febs
        febs = [];
        ix = 1;
        i = 1;
        size = 0;
        while i<=length(evalues)
            multip = repeats(i);
            evalimag = -abs(imag(evalues(i)));
            if -evalimag < precision
                n = 1;
                sigma = -real(evalues(i));
                z = 0;
                ev = real(evalues(i));
                em = multip;
                size = size + 1;
            else
                n = 3;
                size = size + 3;   
                while evalimag / real(evalues(i)) >= cot(pi/n)
                    n = n + 1;
                    size = size + 1;
                    if size > maxSize
                        error ('The represetation is too large (>maxSize). No result returned.');
                    end
                end
                sigma = -(2*real(evalues(i)) + evalimag*(cot(pi/n) - tan(pi/n))) / 2;
                z = (-evalimag*(cot(pi/n) + tan(pi/n)) / (2*sigma))^n;
                ev = []; em = [];
                for k=1:n
                    ev = [ev; complex(-(1-z^(1/n)*cos(2*(k-1)*pi/n))*sigma,z^(1/n)*sin(2*(k-1)*pi/n)*sigma)];
                    em = [em; multip];
                end
            end

            febs(ix).lambda = sigma;
            febs(ix).z = z;
            febs(ix).n = n;
            febs(ix).multip = multip;
            febs(ix).evals = ev;
            febs(ix).emuls = em;

            if -evalimag < precision
                i = i + 1;
            else
                i = i + 2;
            end
            ix = ix + 1;
        end

        % ordering according to the dominant eigenvalue (simple bubble sort)
        N = length(febs);
        for i=1:N-1
            for j=1:N-i
                if max(febs(j).evals) < max(febs(j+1).evals)
                    F = febs(j+1);
                    febs(j+1) = febs(j);
                    febs(j) = F;
                end
            end
        end
    end

    function A = febGenerator (lambda, z, n, multip)
        A = zeros(multip*n);
        for ii=1:multip
            lv = ones(1,n)*lambda;
            A((ii-1)*n+1:ii*n, (ii-1)*n+1:ii*n) = -diag(lv) + diag(lv(1:n-1),1);
            A(ii*n,(ii-1)*n+1) = A(ii*n,(ii-1)*n+1) + z * lambda;
            if ii<multip
                A(ii*n,ii*n+1) = (1-z)*lambda;
            end
        end
    end

    % determine eigenvalues and their multiplicities
    [evalues, repeats] = eigvalc(A);

    % build monocyclic representation
    febs = generateFEBs (A, evalues, repeats, maxSize, precision);

    % assemble generator matrix of the hyper-feb 
    B = [];
    pos = 1;
    for i=1:length(febs)
        Ni = febs(i).n*febs(i).multip;
        B(pos:pos+Ni-1,pos:pos+Ni-1) = febGenerator (febs(i).lambda, febs(i).z, febs(i).n, febs(i).multip);
        if i<length(febs)
            B(pos+Ni-1, pos+Ni) = -sum(B(pos+Ni-1,:));
        end
        pos = pos + Ni;
    end
end
