%  acf = LagCorrelationsFromTrace(trace, K)
%  
%  Returns the lag-k autocorrelation of a trace.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  K : int
%      The number of lags to compute
%  
%  Returns
%  -------
%  acf : column vector of doubles
%      The lag-k autocorrelation function of the trace up to
%      lag K

function acf = LagCorrelationsFromTrace (trace, K)

    if nargin<2
        K = 3;
    end

    m = mean(trace);
    v = var(trace);
    
    acf = zeros(1,K);
    for i=1:K
        acf(i) = (dot(trace(1:end-i),trace(i+1:end)) / (length(trace)-i) - m^2) / v;
    end
end

