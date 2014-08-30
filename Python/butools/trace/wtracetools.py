import numpy as np

def  PdfFromWeightedTrace (trace, weights, intBounds):
    intlens = intBounds[1:] - intBounds[0:-1]
    x = (intBounds[1:] + intBounds[0:-1]) / 2.0
    y = np.empty (len(x))
   
    for i in range(len(x)):
        y[i] = np.sum(weights[np.logical_and(trace>=intBounds[i], trace<intBounds[i+1])])
    
    return (x,y/intlens/np.sum(weights))
  
def  CdfFromWeightedTrace (trace, weights):
    ix = np.argsort(trace)
    return (trace[ix], np.cumsum(weights[ix])/np.sum(weights)) 
  
def MarginalMomentsFromWeightedTrace (trace, weights, K=5):
    return [np.power(trace,k).dot(weights)/np.sum(weights) for k in range(1,K+1)]
  
