from scipy import stats
import numpy as np

def  PdfFromTrace (trace, intBounds):
    hist = stats.histogram2 (trace, intBounds)
    intlens = intBounds[1:] - intBounds[0:-1]
    y = hist[0:-1] / intlens / len(trace)
    x = (intBounds[1:] + intBounds[0:-1]) / 2.0
    return (x,y)
  
def  CdfFromTrace (trace):
    return (np.sort(trace), np.linspace(0.0, 1.0, len(trace)))
  
def IATimesFromCummulative (tr):
    return np.diff(tr)
  
def MarginalMomentsFromTrace (trace, K=5):
    return [np.sum(np.power(trace,k))/len(trace) for k in range(1,K+1)]
  
def LagCorrelationsFromTrace (trace, K=3):
    m = np.mean(trace)
    v = np.var(trace)
    return [ (np.dot(trace[0:-i], trace[i:])/(len(trace)-i) - m*m) / v for i in range(1,K+1)]
  
def LagkJointMomentsFromTrace (trace, K=3, L=1):
    return np.array([ [ np.dot(np.power(trace[:-L],i), np.power(trace[L:],j))/(len(trace)-L) for j in range(0,K+1) ] for i in range(0,K+1)])
