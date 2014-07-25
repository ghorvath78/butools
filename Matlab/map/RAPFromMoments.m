%  [H0, H1] = RAPFromMoments(moms, Nm)
%  
%  Creates a rational arrival process that has the same 
%  marginal and lag-1 joint moments as given (see [1]_).
%  
%  Parameters
%  ----------
%  moms : vector of doubles
%      The list of marginal moments. To obtain a rational 
%      process of order M, 2*M-1 marginal moments are 
%      required.
%  Nm : matrix, shape (M,M)
%      The matrix of lag-1 joint moments. 
%  
%  Returns
%  -------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational process
%  
%  Notes
%  -----
%  There is no guarantee that the returned matrices define
%  a valid stochastic process. The joint densities may be
%  negative.
%  
%  References
%  ----------
%  .. [1] G Horvath, M Telek, "A minimal representation of 
%         Markov arrival processes and a moments matching 
%         method," Performance Evaluation 64:(9-12) pp. 
%         1153-1168. (2007)       

function [H0, H1] = RAPFromMoments (moms, Nm)

    H = MRAPFromMoments (moms, {Nm});
    H0 = H{1};
    H1 = H{2};
end

