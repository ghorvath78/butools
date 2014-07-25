%  [alpha, A] = PH2From3Moments(moms, prec)
%  
%  Returns a PH(2) which has the same 3 moments as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(3)
%    The moments to match
%  prec : double, optional
%    Numerical precision, default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,2)
%    The initial probability vector of the PH(2)
%  A : matrix, shape (2,2)
%    Transient generator matrix of the PH(2)
%  
%  Notes
%  -----
%  Raises an error if the moments are not feasible with
%  a PH(2).
%  
%  References
%  ----------
%  .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
%          discrete and continuous phase-type distributions of
%          second order," in In Proc. of UK Performance 
%          Evaluation Workshop, UKPEW, 2002"

function [alpha, A] = PH2From3Moments (moms, prec)

  if ~exist('prec','var')
    prec = 1e-14;
  end

  m1 = moms(1);
  m2 = moms(2);
  m3 = moms(3);

  % check moment boounds
  m2l = APH2ndMomentLowerBound(m1, 2);
  m3l = APH3rdMomentLowerBound(m1, m2, 2);
  m3u = APH3rdMomentUpperBound(m1, m2, 2);
  
  if m2<m2l
    error('The given second moment is not feasible!');
  end
  if m3<m3l
    error('The given third moment is not feasible (too small)!');
  end
  if m3>m3u
    error('The given third moment is not feasible (too large)!');
  end
    
  % check if we have an exponential distribution
  if abs(m2/m1/m1-2.0) < prec
    alpha = 1;
    A = -1/m1;
    return;
  end
  
  % calculate parameters
  b = 3.0*m1*m2-m3;
  c = 3.0*m2*m2-2.0*m1*m3;
  e = -2.0*m1*m1+m2;
  a = b*b+6.0*c*e;
  if a<0
    a = 0;
  end
  a = sqrt(a);
  if c>0
    lambda1 = (b - a) / c;
    lambda2 = (b + a) / c;
    p = (-b-6.0*m1*e+a) / (b+a);
  elseif c<0
    lambda1 = (b + a) / c;
    lambda2 = (b - a) / c;
    p = (b+6.0*m1*e+a) / (-b+a);
  elseif c==0
    lambda1 = 0;
    lambda2 = 1.0 / m1;
    p = 0;
  end
  
  % return the result
  alpha = [p,1.0-p];
  A = [-lambda1, lambda1; 0,-lambda2];
end
