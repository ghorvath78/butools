%  [alpha, A] = DPH2From3Moments(moms, prec)
%  
%  Returns an order-2 discrete phase-type distribution 
%  which has the same 3 moments as given.
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
%    The initial probability vector of the DPH(2)
%  A : matrix, shape (2,2)
%    Transition probability matrix of the DPH(2)
%  
%  Notes
%  -----
%  Raises an error if the moments are not feasible with
%  a DPH(2).
%  
%  This procedure first calls 'MGFromMoments', then transforms
%  it to DPH(2) by 'CanonicalFromDPH2'.

function [alpha, A] = DPH2From3Moments (moms)

    [beta, B] = MGFromMoments(moms(1:3));
    [alpha,A] = CanonicalFromDPH2(beta,B);    
end
