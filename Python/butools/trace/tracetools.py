from scipy import stats
import numpy as np

def PdfFromTrace (trace, intBounds):
    """
    Returns the empirical density function of a trace.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    intBounds : vector of doubles
        The array of interval boundaries. The pdf is the
        number of samples falling into an interval divided
        by the interval length.
    
    Returns
    -------
    x : vector of doubles
        The center of the intervals (the points where the 
        empirical pdf is calculated)
    y : vector of doubles
        The values of the empirical pdf at the given points
    """
    hist = stats.histogram2 (trace, intBounds)
    intlens = intBounds[1:] - intBounds[0:-1]
    y = hist[0:-1] / intlens / len(trace)
    x = (intBounds[1:] + intBounds[0:-1]) / 2.0
    return (x,y)
  
def CdfFromTrace (trace):
    """
    Returns the empirical distribution function of the trace.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    
    Returns
    -------
    x : vector of doubles
        The points where the empirical cdf is calculated
    y : vector of doubles
        The values of the empirical cdf at the given points
    """
    return (np.sort(trace), np.linspace(0.0, 1.0, len(trace)))
  
def IATimesFromCummulative (tr):
    """
    Returns the vector of inter-arrival times of a trace
    containing cummulative data.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data (cummulative)
    
    Returns
    -------
    iat : vector of doubles
        The inter-arrival times
    """
    return np.diff(tr)
  
def MarginalMomentsFromTrace (trace, K=5):
    """
    Returns the marginal moments of a trace.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The number of moments to compute
    
    Returns
    -------
    moms : vector of doubles
        The (raw) moments of the trace
    """
    return [np.sum(np.power(trace,k))/len(trace) for k in range(1,K+1)]
  
def LagCorrelationsFromTrace (trace, K=3):
    """
    Returns the lag-k autocorrelation of a trace.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The number of lags to compute
    
    Returns
    -------
    acf : column vector of doubles
        The lag-k autocorrelation function of the trace up to
        lag K
    """
    m = np.mean(trace)
    v = np.var(trace)
    return [ (np.dot(trace[0:-i], trace[i:])/(len(trace)-i) - m*m) / v for i in range(1,K+1)]
  
def LagkJointMomentsFromTrace (trace, K=3, L=1):
    """
    Returns the lag-L joint moments of a trace.
    
    It is computed as `Nm_{i,j}=\frac{1}{N-L}\sum_{k=0}^{N-L} x_k^i x_{k+L}^j`.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The joint moments are computed up to order K
    L : int, optional
        The lag-L joint moments are computed.
        The default value is 1.
    
    Returns
    -------
    Nm : matrix, shape (K,K)
        The matrix of joint moments, starting from moment 0
    """
    return np.array([ [ np.dot(np.power(trace[:-L],i), np.power(trace[L:],j))/(len(trace)-L) for j in range(0,K+1) ] for i in range(0,K+1)])
