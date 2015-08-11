%  [x, y] = CdfFromTrace(trace)
%  
%  Returns the empirical distribution function of the trace.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  
%  Returns
%  -------
%  x : vector of doubles
%      The points where the empirical cdf is calculated
%  y : vector of doubles
%      The values of the empirical cdf at the given points

function [x, y] = CdfFromTrace (trace)

    x = reshape(sort(trace), 1, length(trace));
    y = linspace(0, 1, length(trace));
end

