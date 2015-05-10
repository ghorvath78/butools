%  [D0, D1] = DMAP2FromMoments(moms, corr1)
%  
%  Returns a discrete MAP(2) which has the same 3 marginal
%  moments and lag-1 autocorrelation as given.
%  
%  Parameters
%  ----------
%  moms : vector, length(3)
%      First three marginal moments of the inter-arrival times
%  corr1 : double
%      The lag-1 autocorrelation of the inter-arrival times
%  
%  Returns
%  -------
%  D0 : matrix, shape (2,2)
%      The D0 matrix of the discrete MAP(2)
%  D1 : matrix, shape (2,2)
%      The D1 matrix of the discrete MAP(2)
%  
%  Notes
%  -----
%  Raises an exception if the moments are not feasible with
%  a DMAP(2). This procedure calls :func:`butools.dmap.DRAPFromMoments`
%  followed by :func:`butools.dmap.CanonicalFromDMAP2`.
%     

function [D0, D1] = DMAP2FromMoments (moms, corr1)

    Nm = [1, moms(1); moms(1), corr1*(moms(2)-moms(1)^2)+moms(1)^2];
    [H0, H1] = DRAPFromMoments (moms, Nm);

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   
    
    oldCheckInput = BuToolsCheckInput;
    BuToolsCheckInput = false;
    
    [D0, D1] = CanonicalFromDMAP2 (H0, H1);
    
    BuToolsCheckInput = oldCheckInput;
end
