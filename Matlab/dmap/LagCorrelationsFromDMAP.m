%  acf = LagCorrelationsFromDMAP(D0, D1, L, prec)
%  
%  Returns the lag autocorrelations of a discrete Markovian 
%  arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the discrete Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the discrete Markovian arrival process
%  L : double, optional
%      The number of lags to compute. The default value is 1
%  prec : double, optional
%      Numerical precision to check if the input is valid. 
%      The default value is 1e-14
%  
%  Returns
%  -------
%  acf : column vector of doubles, length (L)
%      The lag autocorrelation function up to lag L
%      

function [lagcorrs] = LagCorrelationsFromDMAP (d0, d1, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMAPRepresentation(d0,d1)
        error('LagCorrelationsFromDMAP: input isn''t a valid DMAP representation!');
    end

    if ~exist('L','var')
        L=1;
    end

    lagcorrs = LagCorrelationsFromDRAP(d0, d1, L);
end