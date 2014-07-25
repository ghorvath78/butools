%  iat = IATimesFromCummulative(trace)
%  
%  Returns the vector of inter-arrival times of a trace
%  containing cummulative data.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data (cummulative)
%  
%  Returns
%  -------
%  iat : vector of doubles
%      The inter-arrival times

function iat = IATimesFromCummulative (tr)

    iat = diff (tr);
end

