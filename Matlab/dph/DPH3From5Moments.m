%  [alpha, A] = DPH3From5Moments(moms, prec)
%  
%  Returns an order-3 discrete phase-type distribution 
%  which has the same 5 moments as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(5)
%    The moments to match
%  prec : double, optional
%    Numerical precision, default value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,3)
%    The initial probability vector of the DPH(3)
%  A : matrix, shape (3,3)
%    Transition probability matrix of the DPH(3)
%  
%  Notes
%  -----
%  Raises an error if the moments are not feasible with
%  a DPH(3).
%  
%  This procedure first calls 'MGFromMoments', then transforms
%  it to DPH(3) by 'CanonicalFromDPH3'.

function [alpha, A] = DPH3From5Moments (moms, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end
    
    [beta, B] = MGFromMoments(moms(1:5));
    [alpha,A] = CanonicalFromDPH3(beta,B,prec);    
end
