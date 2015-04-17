%  pi = GM1StationaryDistr (B, R, K)
%  
%  Returns the stationary distribution of the G/M/1 type
%  Markov chain up to a given level K.
%  
%  Parameters
%  ----------
%  A : matrix, shape (N,M*N)
%      Matrix blocks of the G/M/1 type generator in the
%      regular part, from 0 to M-1, concatenated 
%      horizontally.
%  B : matrix, shape (N,M*N)
%      Matrix blocks of the G/M/1 type generator at the
%      boundary, from 0 to M-1, concatenated horizontally.
%  R : matrix, shape (N,N)
%      Matrix R of the G/M/1 type Markov chain
%  K : integer
%      The stationary distribution is returned up to
%      this level.
%  
%  Returns
%  -------
%  pi : matrix, shape (1,(K+1)*N)
%      The stationary probability vector up to level K

function pi = GM1StationaryDistr (B, R, K)

    global BuToolsVerbose;

    pi = GIM1_pi(B,R,'MaxNumComp', K+1, 'Verbose', BuToolsVerbose);
end
