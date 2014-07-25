%  [x, y] = PdfFromWeightedTrace(trace, weights, intBounds)
%  
%  Returns the empirical density function of a trace 
%  consisting of weighted data.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  weights : vector of doubles
%      The weights corresponding to the trace data
%  intBounds : vector of doubles
%      The array of interval boundaries. The pdf is the 
%      number of samples falling into an interval divided by
%      the interval length.
%  
%  Returns
%  -------
%  x : vector of doubles
%      The center of the intervals (the points where the 
%      empirical pdf is calculated)
%  y : vector of doubles
%      The values of the empirical pdf at the given points

function [x, y] = PdfFromWeightedTrace (trace, weights, intBounds)

    intlens = intBounds(2:end) - intBounds(1:end-1);
    x = reshape((intBounds(2:end) + intBounds(1:end-1)) / 2.0, length(intlens), 1);
    y = zeros (length(intlens), 1);
   
    for i=1:length(x)
        y(i) = sum(weights(and(trace>=intBounds(i),trace<intBounds(i+1))));
    end
    y = y ./ reshape(intlens,length(intlens),1) / sum(weights);
end
