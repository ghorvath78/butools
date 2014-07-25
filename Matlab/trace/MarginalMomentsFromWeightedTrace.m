%  moms = MarginalMomentsFromWeightedTrace(trace, weights, K)
%  
%  Returns the marginal moments of a trace consisting of 
%  weighted data.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  weights : vector of doubles
%      The weights corresponding to the trace data
%  K : int
%      The number of moments to compute
%  
%  Returns
%  -------
%  moms : vector of doubles
%      The (raw) moments of the weighted trace

function moms = MarginalMomentsFromWeightedTrace (trace, weights, K)

    if nargin<2
        K = 5;
    end
    
    moms = zeros(1,K);
    for i=1:K
        moms(i) = sum(dot(trace.^i,weights));
    end
    moms = moms / sum(weights);
end
