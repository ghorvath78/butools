import numpy as np
import numpy.matlib as ml
import butools
from butools.dph import MGFromMoments, CheckMGRepresentation
from butools.ph import CanonicalFromPH3
import scipy.linalg as la

def CanonicalFromDPH2 (alpha,A,prec=1e-14):
    """
    Returns the canonical form of an order-2 discrete phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,2)
        Initial vector of the discrete phase-type distribution
    A : matrix, shape (2,2)
        Transition probability matrix of the discrete phase-type
        distribution
    prec : double, optional
      Numerical precision for checking the input, default value
      is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,2)
      The initial probability vector of the canonical form
    B : matrix, shape (2,2)
      Transition probability matrix of the canonical form
    """

    if butools.checkInput and not CheckMGRepresentation (alpha, A, prec):
        raise Exception("CanonicalFromDPH2: Input is not a valid DPH representation!")

    if butools.checkInput and (A.shape[0]!=2 or A.shape[1]!=2):
        raise Exception("CanonicalFromDPH2: Dimension must be 2!")

    ev = la.eigvals(A)
    ix = np.argsort(-np.abs(np.real(ev)))
    lambd = ev[ix]
    e=ml.ones((2,1))
    p1=(alpha*(e-A*e))[0,0]
    if lambd[0]>0 and lambd[1]>0 and lambd[0] != lambd[1]:
        d1=(1-lambd[0])*(1-p1-lambd[1])/(lambd[0]-lambd[1])
        d2=p1-d1
        beta=ml.matrix([[d1*(lambd[0]-lambd[1])/((1-lambd[0])*(1-lambd[1])),(d1+d2)/(1-lambd[1])]])
        B=ml.matrix([[lambd[0],1-lambd[0]],[0,lambd[1]]])
    elif lambd[0]>0 and lambd[0]==lambd[1]:
        d2=p1
        d1=(1-lambd[0])*(1-d2-lambd[0])/lambd[0]
        beta=ml.matrix([[d1*lambd[0]/(1-lambd[0])**2,d2/(1-lambd[0])]])
        B=ml.matrix([[lambd[0],1-lambd[0]],[0,lambd[0]]])
    elif lambd[0]>0:
        d1=(1-lambd[0])*(1-p1-lambd[1])/(lambd[0]-lambd[1])
        d2=p1-d1
        beta=ml.matrix([[(d1*lambd[0]+d2*lambd[1])/((1-lambd[0])*(1-lambd[1])),(d1+d2)*(1-lambd[0]-lambd[1])/((1-lambd[0])*(1-lambd[1]))]])
        B=ml.matrix([[lambd[0]+lambd[1],1-lambd[0]-lambd[1]],[lambd[0]*lambd[1]/(lambd[0]+lambd[1]-1),0]])
    return (np.real(beta),np.real(B))

def CanonicalFromDPH3 (alpha,A,prec=1e-14):
    """
    Returns the canonical form of an order-3 discrete phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,3)
        Initial vector of the discrete phase-type distribution
    A : matrix, shape (3,3)
        Transition probability matrix of the discrete phase-type
        distribution
    prec : double, optional
      Numerical precision for checking the input, default value
      is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,3)
      The initial probability vector of the canonical form
    B : matrix, shape (3,3)
      Transition probability matrix of the canonical form
    """
   
    if butools.checkInput and not CheckMGRepresentation (alpha, A, prec):
        raise Exception("CanonicalFromDPH3: Input is not a valid DPH representation!")

    if butools.checkInput and (A.shape[0]!=3 or A.shape[1]!=3):
        raise Exception("CanonicalFromDPH3: Dimension must be 3!")

    ev = la.eigvals(A)
    ix = np.argsort(-np.abs(np.real(ev)))
    lambd = ev[ix]
    eye = ml.eye(3);

    a0 = -lambd[0]*lambd[1]*lambd[2]
    a1 = lambd[0]*lambd[1]+lambd[0]*lambd[2]+lambd[1]*lambd[2]
    a2 = -lambd[0]-lambd[1]-lambd[2]
    e = ml.matrix([[1,1,1]]).T

    if np.real(lambd[0])>0 and np.real(lambd[1])>=0 and np.real(lambd[2])>=0:
        #PPP case
        alphaout,A2 = CanonicalFromPH3(alpha,A-eye,prec)
        Aout = A2+eye
    elif np.real(lambd[0])>0 and np.real(lambd[1])>=0 and np.real(lambd[2])<0:
        #PPN case
        x1 = lambd[0]
        x2 = lambd[1]+lambd[2]
        x3=lambd[1]*lambd[2]/(lambd[1]+lambd[2]-1)
        Aout=ml.matrix([[x1,1-x1,0],[0,x2,1-x2],[0,x3,0]])
        b3=1/(1-x3)*(e-A*e)
        b2=1/(1-x2)*A*b3
        b1=e-b2-b3
        B=ml.hstack((b1,b2,b3))
        alphaout=alpha*B
    elif np.real(lambd[0])>0 and np.real(lambd[1])<0 and np.real(lambd[2])>=0:
        #PNP case
        x1=-a2
        x2=(a0-a1*a2)/(a2*(1+a2))
        x3=a0*(1+a2)/(a0-a2-a1*a2-a2**2)
        Aout=ml.matrix([[x1,1-x1,0],[x2,0,1-x2],[0,x3,0]])
        b3=1/(1-x3)*(e-A*e)
        b2=1/(1-x2)*A*b3
        b1=e-b2-b3
        if alpha*b1>=0:
            B=ml.hstack((b1,b2,b3))
            alphaout = alpha*B
        else:
            #Set the initial vector first element to 0
            x33=0
            a1=-1
            while x33 <= 1:
                [a1,x1,x2,x3,B]=firstInitElem(x33,lambd,alpha,A)
                if a1 >= 0 and x1 >= 0 and x2 >= 0 and x3 >= 0 and x3+x33 < 1:
                    break
                x33=x33+0.01
            
            if a1 >= 0:
                Aout=ml.matrix([[x1,1-x1,0],[x2,0,1-x2],[0,x3,x33]])
                alphaout=alpha*B
            else:
                #PNP+
                x1=lambd[2]
                x2=lambd[0]+lambd[1]
                x3=lambd[0]*lambd[1]/(lambd[0]+lambd[1]-1)
                Aout=ml.matrix([[x1,0,0],[0,x2,1-x2],[0,x3,0]])

                p1=alpha*(e-A*e)
                p2=alpha*A*(e-A*e)
                d1=(1-lambd[0])*((1-lambd[1])*(1-lambd[2])+(-1+lambd[1]+lambd[2])*p1-p2)/((lambd[0]-lambd[1])*(lambd[0]-lambd[2]))
                d2=(lambd[1]-1)*((1-lambd[0])*(1-lambd[2])+(-1+lambd[0]+lambd[2])*p1-p2)/((lambd[0]-lambd[1])*(lambd[1]-lambd[2]))
                d3=(lambd[2]-1)*((1-lambd[0])*(1-lambd[1])+(-1+lambd[0]+lambd[1])*p1-p2)/((lambd[1]-lambd[2])*(lambd[2]-lambd[0]))
                alphaout=ml.matrix([[d3/(1-lambd[2]),(d1*lambd[0]+d2*lambd[1])/((1-lambd[0])*(1-lambd[1])),(d1+d2)*(1-lambd[0]-lambd[1])/((1-lambd[0])*(1-lambd[1]))]])

                if np.min(alphaout) < 0 or np.min(np.min(Aout)) < 0:
                    raise Exception("DPH3Canonical: Unhandled PNP case!")
    elif np.real(lambd[0])>0 and np.real(lambd[1])<0 and np.real(lambd[2])<0:
        #PNN case
        if np.all(np.isreal(lambd)) or np.abs(lambd[1])**2 <= 2*lambd[0]*(-np.real(lambd[1])):
            x1=-a2
            x2=-a1/(1+a2)
            x3=-a0/(1+a1+a2)
            Aout=ml.matrix([[x1,1-x1,0],[x2,0,1-x2],[x3,0,0]])
            b3=1/(1-x3)*(e-A*e)
            b2=1/(1-x2)*A*b3
            b1=e-b2-b3
            B=ml.hstack((b1,b2,b3))
            alphaout=alpha*B
        else:
           [alphaout,A2]=CanonicalFromPH3(alpha,A-eye,prec)
           Aout=A2+eye
    return (np.real(alphaout),np.real(Aout))

def firstInitElem(m33,sortEigs,alpha,A):

    l1=sortEigs[0]
    l2=sortEigs[1]
    l3=sortEigs[2]
    
    m1=-m33+l1+l2+l3
    m2=-((l2-l3)*(l1**2-l1*l2-l1*l3+l2*l3)*(m33**3-2*m33**2*l1+m33*l1**2-2*m33**2*l2+3*m33*l1*l2-l1**2*l2+m33*l2**2- \
      l1*l2**2-2*m33**2*l3+3*m33*l1*l3-l1**2*l3+3*m33*l2*l3-2*l1*l2*l3-l2**2*l3+m33*l3**2-l1*l3**2-l2*l3**2))/ \
      (2*m33*l1**2*l2+2*m33**2*l1**2*l2-l1**3*l2-3*m33*l1**3*l2+l1**4*l2-2*m33*l1*l2**2-2*m33**2*l1*l2**2+l1**3*l2**2+ \
      l1*l2**3+3*m33*l1*l2**3-l1**2*l2**3-l1*l2**4-2*m33*l1**2*l3-2*m33**2*l1**2*l3+l1**3*l3+3*m33*l1**3*l3-l1**4*l3+ \
      2*m33*l2**2*l3+2*m33**2*l2**2*l3-l2**3*l3-3*m33*l2**3*l3+l2**4*l3+2*m33*l1*l3**2+2*m33**2*l1*l3**2-l1**3*l3**2- \
      2*m33*l2*l3**2-2*m33**2*l2*l3**2+l2**3*l3**2-l1*l3**3-3*m33*l1*l3**3+l1**2*l3**3+l2*l3**3+3*m33*l2*l3**3-l2**2*l3**3+ \
      l1*l3**4-l2*l3**4)
    m3=((l2-l3)*(l1**2-l1*l2-l1*l3+l2*l3)*(m33**3+m33**4-m33**2*l1-2*m33**3*l1+m33**2*l1**2-m33**2*l2-2*m33**3*l2+ \
      m33*l1*l2+3*m33**2*l1*l2-m33*l1**2*l2+m33**2*l2**2-m33*l1*l2**2-m33**2*l3-2*m33**3*l3+m33*l1*l3+3*m33**2*l1*l3- \
      m33*l1**2*l3+m33*l2*l3+3*m33**2*l2*l3-l1*l2*l3-4*m33*l1*l2*l3+l1**2*l2*l3-m33*l2**2*l3+l1*l2**2*l3+m33**2*l3**2- \
      m33*l1*l3**2-m33*l2*l3**2+l1*l2*l3**2))/(-2*m33*l1**2*l2-2*m33**2*l1**2*l2-m33**3*l1**2*l2+l1**3*l2+3*m33*l1**3*l2+ \
      2*m33**2*l1**3*l2-l1**4*l2-m33*l1**4*l2+2*m33*l1*l2**2+2*m33**2*l1*l2**2+m33**3*l1*l2**2-l1**3*l2**2-2*m33*l1**3*l2**2+ \
      l1**4*l2**2-l1*l2**3-3*m33*l1*l2**3-2*m33**2*l1*l2**3+l1**2*l2**3+2*m33*l1**2*l2**3+l1*l2**4+m33*l1*l2**4-l1**2*l2**4+ \
      2*m33*l1**2*l3+2*m33**2*l1**2*l3+m33**3*l1**2*l3-l1**3*l3-3*m33*l1**3*l3-2*m33**2*l1**3*l3+l1**4*l3+m33*l1**4*l3- \
      2*m33*l2**2*l3-2*m33**2*l2**2*l3-m33**3*l2**2*l3+l2**3*l3+3*m33*l2**3*l3+2*m33**2*l2**3*l3-l2**4*l3-m33*l2**4*l3- \
      2*m33*l1*l3**2-2*m33**2*l1*l3**2-m33**3*l1*l3**2+l1**3*l3**2+2*m33*l1**3*l3**2-l1**4*l3**2+2*m33*l2*l3**2+2*m33**2*l2*l3**2+ \
      m33**3*l2*l3**2-l2**3*l3**2-2*m33*l2**3*l3**2+l2**4*l3**2+l1*l3**3+3*m33*l1*l3**3+2*m33**2*l1*l3**3-l1**2*l3**3- \
      2*m33*l1**2*l3**3-l2*l3**3-3*m33*l2*l3**3-2*m33**2*l2*l3**3+l2**2*l3**3+2*m33*l2**2*l3**3-l1*l3**4-m33*l1*l3**4+ \
      l1**2*l3**4+l2*l3**4+m33*l2*l3**4-l2**2*l3**4)
    b3=np.sum(ml.eye(3)-A,1)/(1-m3-m33)
    b2=(-m33*ml.eye(3)+A)*b3/(1-m2)
    b1=(-m3*b3+A*b2)/(1-m1)
    B=ml.hstack((b1,b2,b3))
    a1=alpha*b1
    return (a1,m1,m2,m3,B)
    
def DPH2From3Moments (moms, prec=1e-14):
    """
    Returns an order-2 discrete phase-type distribution 
    which has the same 3 moments as given.
    
    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,2)
      The initial probability vector of the DPH(2)
    A : matrix, shape (2,2)
      Transition probability matrix of the DPH(2)
    
    Notes
    -----
    Raises an error if the moments are not feasible with
    a DPH(2).
    
    This procedure first calls 'MGFromMoments', then transforms
    it to DPH(2) by 'CanonicalFromDPH2'.
    """

    beta, B = MGFromMoments(moms[0:3])
    return CanonicalFromDPH2(beta,B,prec)
    
def DPH3From5Moments (moms, prec=1e-14):
    """
    Returns an order-3 discrete phase-type distribution 
    which has the same 5 moments as given.
    
    Parameters
    ----------
    moms : vector of doubles, length(5)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,3)
      The initial probability vector of the DPH(3)
    A : matrix, shape (3,3)
      Transition probability matrix of the DPH(3)
    
    Notes
    -----
    Raises an error if the moments are not feasible with
    a DPH(3).
    
    This procedure first calls 'MGFromMoments', then transforms
    it to DPH(3) by 'CanonicalFromDPH3'.
    """

    beta, B = MGFromMoments(moms[0:5])
    return CanonicalFromDPH3(beta,B,prec)

