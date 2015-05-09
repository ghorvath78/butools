%  [G0, G1] = CanonicalFromMAP2(D0, D1, prec)
%  
%  Returns the canonical form of an order-2 Markovian
%  arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (2,2)
%      The D0 matrix of the MAP(2)
%  D1 : matrix, shape (2,2)
%      The D1 matrix of the MAP(2)
%  prec : double, optional
%      Numerical precision to check the input, default 
%      value is 1e-14
%  
%  Returns
%  -------
%  G0 : matrix, shape (1,2)
%      The D0 matrix of the canonical MAP(2)
%  G1 : matrix, shape (2,2)
%      The D1 matrix of the canonical MAP(2)
%  
%  Notes
%  -----
%  This procedure calculates 3 marginal moments and the lag-1
%  autocorrelation of the input and calls 'MAP2FromMoments'.

function [G0, G1] = CanonicalFromMAP2 (D0, D1) 

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput
        if size(D0,1)~=2
            error('CanonicalFromMAP2: size is not 2!');
        end
        if ~CheckMAPRepresentation(D0, D1)
	        error('CanonicalFromMAP2: Input isn''t a valid MAP representation!');
        end
    end

    moms = MarginalMomentsFromMAP (D0, D1, 3);
    corr1 = LagCorrelationsFromMAP (D0, D1, 1);
    [G0, G1] = MAP2FromMoments (moms, corr1);
end
