%  x = SamplesFromMAP(D0, D1, K, prec)
%  
%  Generates random samples from a Markovian arrival 
%  process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the Markovian arrival process
%  K : integer
%      The number of samples to generate.
%  prec : double, optional
%      Numerical precision to check if the input Markovian
%      arrival process is valid. The default value is 
%      1e-14.
%  
%  Returns
%  -------
%  x : vector, length(K)
%      The vector of random samples (inter-arrival times).

function x = SamplesFromMAP(D0,D1,k,initial)

    if ~exist('initial','var')
        initial=[];
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1)
        error('SamplesFromMAP: input isn''t a valid MAP representation!');
    end
    
    x = SamplesFromMMAP({D0,D1},k,initial);
end
