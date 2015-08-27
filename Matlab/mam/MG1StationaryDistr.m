%  pi = MG1StationaryDistr (A, B, G, K)
%  
%  Returns the stationary distribution of the M/G/1 type
%  Markov chain up to a given level K.
%  
%  Parameters
%  ----------
%  A : length(M) list of matrices of shape (N,N)
%      Matrix blocks of the M/G/1 type generator in the 
%      regular part, from 0 to M-1.
%  B : length(M) list of matrices of shape (N,N)
%      Matrix blocks of the M/G/1 type generator at the
%      boundary, from 0 to M-1.
%  G : matrix, shape (N,N)
%      Matrix G of the M/G/1 type Markov chain
%  K : integer
%      The stationary distribution is returned up to
%      this level.
%  
%  Returns
%  -------
%  pi : array, shape (1,(K+1)*N)
%      The stationary probability vector up to level K

function pi = MG1StationaryDistr (A, B, G, K)

    global BuToolsVerbose;

    N = size(A{1},2);
    Am = zeros(N, N*length(A));
    for i=1:length(A)
        Am(:,(i-1)*N+1:i*N) = A{i};
    end

    Bm = zeros(N, N*length(B));
    for i=1:length(B)
        Bm(:,(i-1)*N+1:i*N) = B{i};
    end

    pi = MG1_pi(Bm,Am,G,'MaxNumComp', K+1, 'Verbose', double(BuToolsVerbose));
end
