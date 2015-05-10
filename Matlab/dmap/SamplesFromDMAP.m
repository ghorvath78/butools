%  x = SamplesFromDMAP(D0, D1, K, prec)
%  
%  Generates random samples from a discrete Markovian 
%  arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the discrete MAP.
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the discrete MAP.
%  K : integer
%      The number of samples to generate.
%  prec : double, optional
%      Numerical precision to check if the input DMAP is
%      valid. The default value is 1e-14.
%  
%  Returns
%  -------
%  x : vector, length(K)
%      The vector of random samples (inter-arrival times).

function x = SamplesFromDMAP(D0,D1,k,initial,prec)

    if ~exist('initial','var')
        initial=[];
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMAPRepresentation(D0,D1)
        error('SamplesFromDMAP: input isn''t a valid DMAP representation!');
    end
    
    x = SamplesFromDMMAP({D0,D1},k,initial);
end
