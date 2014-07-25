%  m3 = APH3rdMomentUpperBound(m1, m2, n)
%  
%  Returns the upper bound of the third moment of acyclic
%  phase-type (APH) distributions of order n.
%  
%  Parameters
%  ----------
%  m1 : double
%      The first moment
%  m2 : double
%      The second moment
%  n : int
%      Number of states
%  
%  Returns
%  -------
%  m3 : double
%      The highest third moment an order-n APH can have with
%      the given first and second moment.
%  
%  References
%  ----------
%  .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
%         moments with minimal acyclic phase type 
%         distributions," Stochastic models, pp. 303-326, 2005.

function m3 = APH3rdMomentUpperBound (m1, m2, n)

    n2 = m2 / m1 / m1;
    if n2<(n+1.0)/n
        m3 = -inf;
    elseif n2<=n/(n-1.0)
        m3 = m1 * m2 * (2.0*(n-2.0)*(n*n2-n-1.0)*sqrt(1.0+(n*(n2-2.0))/(n-1.0)) + (n+2.0)*(3.0*n*n2-2.0*n-2.0)) / (n*n*n2);
    else
        m3 = inf;
    end
end
