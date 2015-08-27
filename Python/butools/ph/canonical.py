import numpy as np
import numpy.matlib as ml
import math
import cmath
import butools
from butools.ph import CheckMERepresentation, MomentsFromPH
from butools.moments import ReducedMomsFromMoms, NormMomsFromMoms

def APH2ndMomentLowerBound (m1, n):
    """
    Returns the lower bound of the second moment of acyclic 
    phase-type (APH) distributions of order n.
    
    Parameters
    ----------
    m1 : double
        The first moment
    n : int
        Number of states
    
    Returns
    -------
    m2 : double
        The lowest second moment an order-n APH can have with
        the given first moment.
    
    References
    ----------
    .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
            discrete and continuous phase-type distributions of
            second order," in In Proc. of UK Performance 
            Evaluation Workshop, UKPEW, 2002"
    """

    return float(m1)*m1*(n+1) / n
    
def APH3rdMomentLowerBound (m1, m2, n):
    """
    Returns the lower bound of the third moment of acyclic
    phase-type (APH) distributions of order n.
    
    Parameters
    ----------
    m1 : double
        The first moment
    m2 : double
        The second moment
    n : int
        Number of states
    
    Returns
    -------
    m3 : double
        The lowest third moment an order-n APH can have with
        the given first and second moment.
    
    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.
    """

    n2 = m2 / m1 / m1
    if n2<(n+1.0)/n:
        return np.inf
    elif n2<(n+4.0)/(n+1.0):
        p = ((n+1.0)*(n2-2.0)) / (3.0*n2*(n-1.0)) * ((-2.0*math.sqrt(n+1.0)) / cmath.sqrt(-3.0*n*n2+4.0*n+4.0) -1.0)
        a = (n2-2.0) / (p*(1.0-n2) + cmath.sqrt(p*p+p*n*(n2-2.0)/(n-1.0)))
        l = ((3.0+a)*(n-1.0)+2.0*a) / ((n-1.0)*(1.0+a*p)) - (2.0*a*(n+1.0)) / (2.0*(n-1.0)+a*p*(n*a+2.0*n-2.0))
        return l.real * m1 * m2
    else:
        return (n+1.0)/n * n2 * m1 * m2
    
def APH3rdMomentUpperBound (m1, m2, n):
    """
    Returns the upper bound of the third moment of acyclic
    phase-type (APH) distributions of order n.
    
    Parameters
    ----------
    m1 : double
        The first moment
    m2 : double
        The second moment
    n : int
        Number of states
    
    Returns
    -------
    m3 : double
        The highest third moment an order-n APH can have with
        the given first and second moment.
    
    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.
    """

    n2 = m2 / m1 / m1
    if n2<(n+1.0)/n:
        return -np.inf
    elif n2<=n/(n-1.0):
        return m1 * m2 * (2.0*(n-2.0)*(n*n2-n-1.0)*math.sqrt(1.0+(n*(n2-2.0))/(n-1.0)) + (n+2.0)*(3.0*n*n2-2.0*n-2.0)) / (n*n*n2)
    else:
        return np.inf

def APHFrom2Moments (moms):
    """
    Returns an acyclic PH which has the same 2 moments as
    given. If detects the order and the structure 
    automatically to match the given moments.
    
    Parameters
    ----------
    moms : vector of doubles, length(2)
      The moments to match
    maxSize : int, optional
      The maximal size of the resulting APH. The default value
      is 100.
    
    Returns
    -------
    alpha : vector, shape (1,M)
      The initial probability vector of the APH
    A : matrix, shape (M,M)
      Transient generator matrix of the APH
    
    Raises an error if the moments are not feasible with an
    APH of size "maxSize".
    """
    
    m1, m2 = moms
    cv2 = m2/m1/m1 - 1.0
    lambd = 1.0 / m1
    N = max(int(math.ceil(1.0/cv2)), 2)
    p = 1.0 / (cv2 + 1.0 + (cv2-1.0)/(N-1))
    A = -lambd * p * N * ml.eye(N)
    for i in range(0,N-1):
        A[i,i+1] = -A[i,i]
    A[-1,-1] = -lambd * N
    alpha = ml.zeros((1,N))
    alpha[0,0] = p
    alpha[0,N-1] = 1.0-p
    return (alpha,A)

    
def PH2From3Moments (moms, prec=1e-14):
  """
  Returns a PH(2) which has the same 3 moments as given.
  
  Parameters
  ----------
  moms : vector of doubles, length(3)
    The moments to match
  prec : double, optional
    Numerical precision, default value is 1e-14
  
  Returns
  -------
  alpha : matrix, shape (1,2)
    The initial probability vector of the PH(2)
  A : matrix, shape (2,2)
    Transient generator matrix of the PH(2)
  
  Notes
  -----
  Raises an error if the moments are not feasible with
  a PH(2).
  
  References
  ----------
  .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
          discrete and continuous phase-type distributions of
          second order," in In Proc. of UK Performance 
          Evaluation Workshop, UKPEW, 2002"
  """

  m1, m2, m3 = moms

  # check moment boounds
  m2l = APH2ndMomentLowerBound(m1, 2)  
  m3l = APH3rdMomentLowerBound(m1, m2, 2)  
  m3u = APH3rdMomentUpperBound(m1, m2, 2)  
  
  if m2<m2l:
    raise Exception("The given second moment is not feasible!")    
  if m3<m3l:
    raise Exception("The given third moment is not feasible (too small)!")
  if m3>m3u:
    raise Exception("The given third moment is not feasible (too large)!")
    
  # check if we have an exponential distribution
  if abs(m2/m1/m1-2.0) < prec:
    return (np.matrix([1]), np.matrix([[-1/m1]]))
  
  # calculate parameters
  b = 3.0*m1*m2-m3
  c = 3.0*m2*m2-2.0*m1*m3
  e = -2.0*m1*m1+m2
  a = b*b+6.0*c*e
  if a<0:
    a = 0
  a = math.sqrt(a)
  if c>0:
    lambda1 = (b - a) / c
    lambda2 = (b + a) / c
    p = (-b-6.0*m1*e+a) / (b+a)
  elif c<0:
    lambda1 = (b + a) / c
    lambda2 = (b - a) / c
    p = (b+6.0*m1*e+a) / (-b+a)
  elif c==0:
    lambda1 = 0
    lambda2 = 1.0 / m1
    p = 0
  
  # return the result
  return (np.matrix([p,1.0-p]), np.matrix([[-lambda1, lambda1], [0,-lambda2]]))
  
def APHFrom3Moments (moms, maxSize=100, prec=1e-14):
  """
  Returns an acyclic PH which has the same 3 moments as
  given. If detects the order and the structure 
  automatically to match the given moments.
  
  Parameters
  ----------
  moms : vector of doubles, length(3)
    The moments to match
  maxSize : int, optional
    The maximal size of the resulting APH. The default value
    is 100.
  
  Returns
  -------
  alpha : vector, shape (1,M)
    The initial probability vector of the APH
  A : matrix, shape (M,M)
    Transient generator matrix of the APH
  
  Raises an error if the moments are not feasible with an
  APH of size "maxSize".
  
  References
  ----------
  .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
         moments with minimal acyclic phase type 
         distributions," Stochastic models, pp. 303-326, 2005.
  """

  m1, m2, m3 = moms

  # detect number of phases needed
  n = 2
  while n<maxSize and (APH2ndMomentLowerBound(m1, n) > m2 or APH3rdMomentLowerBound(m1, m2, n) >= m3 or APH3rdMomentUpperBound(m1, m2, n) <= m3):
    n = n + 1

  # if PH is too large, adjust moment to bounds
  if APH2ndMomentLowerBound(m1, n) > m2:
    m2 = APH2ndMomentLowerBound(m1, n)

  if APH3rdMomentLowerBound(m1, m2, n) > m3:
    m3 = APH3rdMomentLowerBound(m1, m2, n)

  if APH3rdMomentUpperBound(m1, m2, n) < m3:
    m3 = APH3rdMomentUpperBound(m1, m2, n)

  # compute normalized moments
  n1,n2,n3 = NormMomsFromMoms ([m1, m2, m3])

  if n2>2.0 or n3 < 2.0*n2 - 1.0:
    b = (2.0*(4.0-n*(3.0*n2-4.0)) / (n2*(4.0+n-n*n3) + math.sqrt(n*n2)*math.sqrt(12.0*n2*n2*(n+1.0)+16.0*n3*(n+1.0)+n2*(n*(n3-15.0)*(n3+1.0)-8.0*(n3+3.0))))).real
    a = (b*n2-2.0)*(n-1.0)*b / (b-1.0) / n
    p = (b-1.0) / a
    lamb = (p*a+1.0) / n1
    mu = (n-1.0)*lamb / a
    # construct representation
    alpha = ml.zeros((1,n))
    alpha[0,0] = p
    alpha[0,n-1] = 1.0-p
    A = ml.zeros((n,n))
    A[n-1,n-1] = -lamb;
    for i in range(n-1):
      A[i,i] = -mu
      A[i,i+1] = mu
    return (alpha,A)
  else:
    c4 = n2*(3.0*n2-2.0*n3)*(n-1.0)*(n-1.0)
    c3 = 2.0*n2*(n3-3.0)*(n-1.0)*(n-1.0)
    c2 = 6.0*(n-1.0)*(n-n2)
    c1 = 4.0*n*(2.0-n)
    c0 = n*(n-2.0)
    fs = np.roots([c4, c3, c2, c1, c0])
    for f in fs:
      if abs((n-1)*(n2*f*f*2-2*f+2)-n)<1e-14:
        continue
      a = 2.0*(f-1.0)*(n-1.0) / ((n-1.0)*(n2*f*f-2.0*f+2.0)-n)
      p = (f-1.0)*a
      lamb = (a+p) / n1
      mu = (n-1.0) / (n1 - p/lamb)
      if np.isreal(p) and np.isreal(lamb) and np.isreal(mu) and p>=0 and p<=1 and lamb>0 and mu>0:
        alpha = ml.zeros((1,n))
        alpha[0,0] = p.real
        alpha[0,1] = 1.0-p.real
        A = ml.zeros((n,n));
        A[0,0] = -lamb.real
        A[0,1] = lamb.real
        for i in range(1,n):
          A[i,i] = -mu.real
          if i<n-1:
            A[i,i+1] = mu.real
        return (alpha,A)
  raise Exception("No APH found for the given 3 moments!")
    
def PH3From5Moments (moms, prec=1e-10):
    """
    Returns a PH(3) which has the same 5 moments as given.
    
    Parameters
    ----------
    moms : vector of doubles, length(5)
      The moments to match
    
    Returns
    -------
    alpha : vector, shape (1,3)
        The initial probability vector of the PH(3)
    A : matrix, shape (3,3)
        Transient generator matrix of the PH(3)
    
    Notes
    -----
    Raises an error if the moments are not feasible with
    a PH(3). Also note that the numerical behavior of the 
    procedure can be poor if the moments are close to the 
    boundary of the feasible region.
    
    References
    ----------
    .. [1] G. Horvath and M. Telek, "On the canonical 
           representation of phase type distributions," 
           Performance Evaluation, vol. 66, no. 8, pp. 
           396 - 409, 2009.
    """
    
    m1, m2, m3, m4, m5 = moms
    
    #convert the moments to reduced moments
    moms = ReducedMomsFromMoms([m1,m2,m3,m4,m5])
    for i in range(5):
      moms[i] /= m1**(i+1)
    #solve linear system of equations for a0 a1 a2
    M = np.matrix ([[moms[2], -moms[1], moms[0]],[moms[3], -moms[2], moms[1]],[moms[4], -moms[3], moms[2]]])
    a = np.linalg.solve(M, [1, moms[0], moms[1]])

    discr = a[2]*a[2]-3.0*a[1]
    if discr < 0:
      raise Exception("Invalid characteristic polynomial!")

    gu = (a[2] + 2.0*math.sqrt(discr)) / 3.0
    g0 = (a[2] + math.sqrt(discr)) / 3.0
   
    rts = np.roots(np.hstack((1,a[::-1])))
    ix = np.argsort(np.real(rts))   
    lamb = (-rts)[ix]
    
    d1 = a[1] - a[2] - a[0] * moms[1];
    d2 = a[0] - a[1] - a[2] * d1;
    d3 = -a[0] - a[1]*d1 - a[2]*d2;

    if d1>1e-10 or (abs(d1)<prec and d2>0):
      raise Exception("Negative density around 0!")

    if lamb[2]<0:
      raise Exception("Invalid eigenvalues!")

    if lamb[0].imag < prec:
        gl = lamb[0].real
    else:
        gl = g0

    if gl > gu+prec:
      raise Exception("Invalid eigenvalues (gl>gu detected)!")
    if gl > gu:
        gl = gu

    if abs(d1) < prec:
        g2 = 0
    else:
        g2 = -d2 / d1

    if g2>gu+prec:
      raise Exception("alpha_2 is negative!")
    if g2 > gu:
        g2 = gu;

    x1 = max(g2, gl)

    if np.isreal(lamb[0]) and g2<gl:
      x13 = 0
    else:
      x13 = x1 - a[0] / (x1*x1 - a[2]*x1 + a[1])

    bels = (a[2]-x1)**2 - 4.0*(x1*x1-a[2]*x1+a[1])
    if bels<0 and bels>-prec:
        bels = 0

    x2 = (a[2] - x1 + math.sqrt(bels)) / 2.0
    x3 = (a[2] - x1 - math.sqrt(bels)) / 2.0
    p1 = d1 / (x13 - x1)
    p2 = (x1*d1 + d2) / (x13-x1) / x2
    p3 = (x1*x2*d1 + x2*d2 + x1*d2 + d3) / (x13-x1) / x2 / x3

    T = np.matrix([[-x1, 0, x13],[x2, -x2, 0],[0, x3, -x3]]) / m1
    alpha = np.matrix([p1,p2,p3])

    if x13<-prec or x13>x1:
      raise Exception("Invalid genrator!")

    if np.min(np.real(alpha))<-prec:
      raise Exception("Initial vector has negative entries!")
    
    if np.max(abs(np.imag(alpha)))>prec:
      raise Exception("Inital vector has complex entries!")

    if np.max(np.real(alpha))>1+prec:
      raise Exception("Initial vector has entries that are greater than 1!")

    return (alpha, T)

def CanonicalFromPH2 (alpha, A, prec=1e-14):
    """
    Returns the canonical form of an order-2 phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,2)
      Initial vector of the phase-type distribution
    A : matrix, shape (2,2)
      Transient generator of the phase-type distribution
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,2)
      The initial probability vector of the canonical form
    B : matrix, shape (2,2)
      Transient generator matrix of the canonical form
    
    Notes
    -----
    This procedure calculates 3 moments of the input and
    calls 'PH2From3Moments'.
    """

    if butools.checkInput and not CheckMERepresentation(alpha,A):
        raise Exception("CanonicalFromPH2: Input isn''t a valid ME distribution!")

    if A.shape[0]!=2:
        raise Exception("CanonicalFromPH2: Dimension is not 2!")

    return PH2From3Moments (MomentsFromPH(alpha, A, 3), prec)

def CanonicalFromPH3 (alpha, A, prec=1e-14):
    """
    Returns the canonical form of an order-3 phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,3)
        Initial vector of the phase-type distribution
    A : matrix, shape (3,3)
        Transient generator of the phase-type distribution
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,3)
      The initial probability vector of the canonical form
    B : matrix, shape (3,3)
      Transient generator matrix of the canonical form
    
    Notes
    -----
    This procedure calculates 5 moments of the input and
    calls 'PH3From5Moments'.
    """

    if butools.checkInput and not CheckMERepresentation(alpha,A):
        raise Exception("CanonicalFromPH3: Input isn''t a valid ME distribution!")

    if A.shape[0]!=3:
        raise Exception("CanonicalFromPH3: Dimension is not 3!")
    
    return PH3From5Moments (MomentsFromPH(alpha, A, 5), prec)

