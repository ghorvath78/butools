%  moms = MarginalMomentsFromTrace(trace, K)
%  
%  Returns the marginal moments of a trace.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  K : int
%      The number of moments to compute
%  
%  Returns
%  -------
%  moms : vector of doubles
%      The (raw) moments of the trace

function moms = MarginalMomentsFromTrace (trace, K)

    if nargin<2
        K = 5;
    end
    
    moms = zeros(1,K);
    for i=1:K
        moms(i) = sum(trace.^i);
    end
    moms = moms / length(trace);
end
