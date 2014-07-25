%  m3 = APH3rdMomentLowerBound(m1, m2, n)
%  
%  Returns the lower bound of the third moment of acyclic
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
%      The lowest third moment an order-n APH can have with
%      the given first and second moment.
%  
%  References
%  ----------
%  .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
%         moments with minimal acyclic phase type 
%         distributions," Stochastic models, pp. 303-326, 2005.

function m3 = APH3rdMomentLowerBound (m1, m2, n)

    n2 = m2 / m1 / m1;
    if n2<(n+1.0)/n
        m3 = inf;
    elseif n2<(n+4.0)/(n+1.0)
        p = ((n+1.0)*(n2-2.0)) / (3.0*n2*(n-1.0)) * ((-2.0*sqrt(n+1.0)) / sqrt(-3.0*n*n2+4.0*n+4.0) -1.0);
        a = (n2-2.0) / (p*(1.0-n2) + sqrt(p*p+p*n*(n2-2.0)/(n-1.0)));
        l = ((3.0+a)*(n-1.0)+2.0*a) / ((n-1.0)*(1.0+a*p)) - (2.0*a*(n+1.0)) / (2.0*(n-1.0)+a*p*(n*a+2.0*n-2.0));
        m3 = real(l) * m1 * m2;
    else
        m3 = (n+1.0)/n * n2 * m1 * m2;
    end
end

