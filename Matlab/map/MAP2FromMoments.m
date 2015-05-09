%  [D0, D1] = MAP2FromMoments(moms, corr1)
%  
%  Returns a MAP(2) which has the same 3 marginal moments 
%  and lag-1 autocorrelation as given.
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
%      The D0 matrix of the MAP(2)
%  D1 : matrix, shape (2,2)
%      The D1 matrix of the MAP(2)
%  
%  Raises an exception if the moments are not feasible with
%  a MAP(2).
%  
%  Notes
%  -----
%  The result is always a valid MAP(2) as long as the input
%  moments can be realized by a PH(2) (can be checked with 
%  :func:`butools.ph.APH2ndMomentLowerBound`, 
%  :func:`butools.ph.APH3rdMomentLowerBound`, 
%  :func:`butools.ph.APH3rdMomentUpperBound` with n=2) and the 
%  correlation falls between the feasible lower and upper 
%  bound (check by :func:`MAP2CorrelationBounds`).
%  
%  References
%  ----------
%  .. [1] L Bodrog, A Heindl, G Horvath, M Telek, "A Markovian
%         Canonical Form of Second-Order Matrix-Exponential 
%         Processes," EUROPEAN JOURNAL OF OPERATIONAL RESEARCH
%         190:(2) pp. 459-477. (2008)
%      

function [D0, D1] = MAP2FromMoments (moms, corr1)

    m1 = moms(1);
    m2 = moms(2);
    m3 = moms(3);

    % If we have an exponential distribution, we do not allow correlation
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-12;
    end   
    if abs(m2-2.0*m1*m1) < BuToolsCheckPrecision && abs(corr1) > BuToolsCheckPrecision
        error('We do not allow correlation in case of exponentially distributed marginal');
    end

    % Perform PH fitting  
    [tau, T] =  PH2From3Moments (moms);
    l1 = -T(1,1);
    l2 = -T(2,2);
    p = tau(1);
    alpha = l1/l2;

    % Check the feasibility of the correlation parameter
    [corrl, corru] = MAP2CorrelationBounds (moms);
    if corr1 < corrl
        error('The correlation parameter is too small!');
    end
    if corr1 > corru
        error('The correlation parameter is too large!');
    end

    gamma = corr1 * (m2-m1*m1) / (m2/2.0-m1*m1);

    % Perform MAP fitting
    if gamma > 0
        a = (1.0+alpha*gamma-p*(1.0-gamma)-sqrt((1.0+alpha*gamma-p*(1.0-gamma))^2-4.0*alpha*gamma)) / (2.0*alpha);
        b = (1.0+alpha*gamma-p*(1.0-gamma)+sqrt((1.0+alpha*gamma-p*(1.0-gamma))^2-4.0*alpha*gamma)) / 2.0;
        D0 = [-l1, (1.0-a)*l1; 0, -l2];
        D1 = [a*l1, 0; (1.0-b)*l2, b*l2];
    elseif gamma < 0
        a = gamma / (alpha*gamma-p*(1.0-gamma));
        b = p*(1.0-gamma)-alpha*gamma;
        D0 = [-l1, (1.0-a)*l1; 0, -l2];
        D1 = [0, a*l1; b*l2, (1.0-b)*l2];
    else
        D0 = [-l1, l1; 0, -l2];
        D1 = [0, 0; p*l2, (1.0-p)*l2];
    end
end
