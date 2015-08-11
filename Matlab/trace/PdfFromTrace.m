%  [x, y] = PdfFromTrace(trace, intBounds)
%  
%  Returns the empirical density function of a trace.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  intBounds : vector of doubles
%      The array of interval boundaries. The pdf is the
%      number of samples falling into an interval divided
%      by the interval length.
%  
%  Returns
%  -------
%  x : vector of doubles
%      The center of the intervals (the points where the 
%      empirical pdf is calculated)
%  y : vector of doubles
%      The values of the empirical pdf at the given points

function [x, y] = PdfFromTrace (trace, intBounds)

    intBounds = reshape(intBounds,length(intBounds),1);
    hist = histc (trace, intBounds);
    intlens = intBounds(2:end) - intBounds(1:end-1);
    y = hist(1:end-1) ./ intlens / length(trace);
    x = (intBounds(2:end) + intBounds(1:end-1)) / 2.0;
    y = reshape(y, 1, length(y));
    x = reshape(x, 1, length(x));   
end
