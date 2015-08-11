%  acf = LagCorrelationsFromRAP(H0, H1, L, prec)
%  
%  Returns the lag autocorrelations of a rational arrival
%  process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
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

function acf = LagCorrelationsFromRAP (H0, H1, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckRAPRepresentation(H0,H1)
        error('LagCorrelationsFromRAP: input isn''t a valid MAP representation!');
    end

    if ~exist('L','var')
        L=1;
    end

    H0i = inv(-H0);
    P = H0i*H1;
    pi = DRPSolve(P);
    moms = MomentsFromME(pi, H0, 2);
    pi = pi * H0i * P;

    acf = zeros(1,L);
    for i=1:L
        acf(i) = (sum(pi*H0i) - moms(1)^2) / (moms(2) - moms(1)^2);
        pi = pi * P;
    end
end

