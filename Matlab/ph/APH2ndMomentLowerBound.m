%  m2 = APH2ndMomentLowerBound(m1, n)
%  
%  Returns the lower bound of the second moment of acyclic 
%  phase-type (APH) distributions of order n.
%  
%  Parameters
%  ----------
%  m1 : double
%      The first moment
%  n : int
%      Number of states
%  
%  Returns
%  -------
%  m2 : double
%      The lowest second moment an order-n APH can have with
%      the given first moment.
%  
%  References
%  ----------
%  .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
%          discrete and continuous phase-type distributions of
%          second order," in In Proc. of UK Performance 
%          Evaluation Workshop, UKPEW, 2002"

function m2 = APH2ndMomentLowerBound (m1, n)

    m2 = m1*m1*(n+1) / n;
end

