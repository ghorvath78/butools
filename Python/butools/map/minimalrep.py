# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:01:03 2013

@author: gabor
"""

import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.map.basemap import MarginalDistributionFromMRAP
from butools.reptrans import MStaircase
from butools.map import CheckRAPRepresentation, CheckMRAPRepresentation

def MinimalRepFromMRAP (H, how="obscont", precision=1e-12):
    """
    Returns the minimal representation of a marked rational
    arrival process.
    
    Parameters
    ----------
    H : list of matrices of shape (M,M)
        The list of H0, H1, ..., HK matrices of the marked
        rational arrival process
    how : {"obs", "cont", "obscont"}, optional        
        Determines how the representation is minimized. 
        "cont" means controllability, "obs" means 
        observability, "obscont" means that the rational arrival
        process is minimized in both respects. Default value 
        is "obscont".
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.
    
    Returns
    -------
    D : list of matrices of shape (M,M)
        The D0, D1, ..., DK matrices of the minimal 
        representation
    
    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.
    """

    if butools.checkInput and not CheckMRAPRepresentation (H):
        raise Exception("MinimalRepFromMRAP: Input is not a valid MRAP representation!")    

    if how=="cont":
        B, n = MStaircase (H, ml.ones((H[0].shape[0],1)), precision)
        return [(la.inv(B)*Hk*B)[0:n,0:n] for Hk in H]
    elif how=="obs":
        alpha, A = MarginalDistributionFromMRAP (H)
        G = [Hk.T for Hk in H]
        B, n = MStaircase (G, alpha.T, precision)
        return [(la.inv(B)*Hk*B)[0:n,0:n] for Hk in H]
    elif how=="obscont":
        D = MinimalRepFromMRAP(H, "cont", precision)
        return MinimalRepFromMRAP(D, "obs", precision)
   
def MinimalRepFromRAP (H0, H1, how="obscont", precision=1e-12):
    """
    Returns the minimal representation of a rational arrival 
    process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    how : {"obs", "cont", "obscont"}, optional      
        Determines how the representation is minimized. "cont" 
        means controllability, "obs" means observability, 
        "obscont" means that the rational arrival process is 
        minimized in both respects. The default value is 
        "obscont"
    precision : double, optional
       Precision used by the Staircase algorithm. The default 
       value is 1e-12.
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the minimal representation
    D1 : matrix, shape (M,M)
        The D1 matrix of the minimal representation
    
    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.
    """

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("MinimalRepFromRAP: Input is not a valid RAP representation!")    

    return MinimalRepFromMRAP ([H0,H1], how, precision)
