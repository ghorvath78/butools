%  pi = MG1StationaryDistr (A, B, G, K)
%  
%  Returns the stationary distribution of the M/G/1 type
%  Markov chain up to a given level K.
%  
%  Parameters
%  ----------
%  A : matrix, shape (N,M*N)
%      Matrix blocks of the M/G/1 type generator in the
%      regular part, from 0 to M-1, concatenated 
%      horizontally.
%  B : matrix, shape (N,M*N)
%      Matrix blocks of the M/G/1 type generator at the
%      boundary, from 0 to M-1, concatenated horizontally.
%  G : matrix, shape (N,N)
%      Matrix G of the M/G/1 type Markov chain
%  K : integer
%      The stationary distribution is returned up to
%      this level.
%  
%  Returns
%  -------
%  pi : matrix, shape (1,(K+1)*N)
%      The stationary probability vector up to level K

function pi = MG1StationaryDistr (A, B, G, K)

    global BuToolsVerbose;

    pi = MG1_pi(B,A,G,'MaxNumComp', K+1, 'Verbose', BuToolsVerbose);
end
