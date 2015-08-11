%  [x, y] = CdfFromWeightedTrace(trace, weights)
%  
%  Returns the empirical distribution function of a trace
%  consisting of weighted data.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  weights : vector of doubles
%      The weights corresponding to the trace data
%  
%  Returns
%  -------
%  x : vector of doubles
%      The points where the empirical cdf is calculated
%  y : vector of doubles
%      The values of the empirical cdf at the given points

function [x, y] = CdfFromWeightedTrace (trace, weights)

    [x,ix] = sort (trace);
    y = cumsum(weights(ix))/sum(weights);
    x = reshape(x, 1, length(x));
    y = reshape(y, 1, length(y));
end
