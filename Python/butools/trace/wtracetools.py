import numpy as np

def PdfFromWeightedTrace (trace, weights, intBounds):
    """
    Returns the empirical density function of a trace 
    consisting of weighted data.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    intBounds : vector of doubles
        The array of interval boundaries. The pdf is the 
        number of samples falling into an interval divided by
        the interval length.
    
    Returns
    -------
    x : vector of doubles
        The center of the intervals (the points where the 
        empirical pdf is calculated)
    y : vector of doubles
        The values of the empirical pdf at the given points
    """
    intlens = intBounds[1:] - intBounds[0:-1]
    x = (intBounds[1:] + intBounds[0:-1]) / 2.0
    y = np.empty (len(x))
   
    for i in range(len(x)):
        y[i] = np.sum(np.array(weights)[np.logical_and(np.array(trace)>=intBounds[i], np.array(trace)<intBounds[i+1])])
    
    return (x,y/intlens/np.sum(weights))
  
def CdfFromWeightedTrace (trace, weights):
    """
    Returns the empirical distribution function of a trace
    consisting of weighted data.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    
    Returns
    -------
    x : vector of doubles
        The points where the empirical cdf is calculated
    y : vector of doubles
        The values of the empirical cdf at the given points
    """
    ix = np.argsort(trace)
    return (np.array(trace)[ix], np.cumsum(np.array(weights)[ix])/np.sum(weights)) 
  
def MarginalMomentsFromWeightedTrace (trace, weights, K=5):
    """
    Returns the marginal moments of a trace consisting of 
    weighted data.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    K : int
        The number of moments to compute
    
    Returns
    -------
    moms : vector of doubles
        The (raw) moments of the weighted trace
    """
    return [np.power(trace,k).dot(weights)/np.sum(weights) for k in range(1,K+1)]
  
