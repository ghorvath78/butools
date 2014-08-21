import numpy as np
import scipy.linalg as la

def NormMomsFromMoms (m):
    """
    Returns the normalized moments given the raw moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The normalized moments are: `\displaystyle n_i=\frac{m_i}{m_{i-1} m_1}`
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    nm : vector of doubles
        The list of normalized moments
    """
    return [float(m[i])/m[i-1]/m[0] if i>0 else m[0] for i in range(len(m))]

def MomsFromNormMoms (nm):
    """
    Returns the raw moments given the normalized moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The normalized moments are: `\displaystyle n_i=\frac{m_i}{m_{i-1} m_1}`
       
    Parameters
    ----------
    nm : vector of doubles
        The list of normalized moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments
    """
    m = nm[:]
    for i in range(1,len(nm)):
        m[i] *= m[0] * m[i-1]
    return m

def ReducedMomsFromMoms (m):
    """
    Returns the reduced moments given the raw moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The reduced moments are: `\displaystyle r_i=\frac{m_i}{i!}`
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    rm : vector of doubles
        The list of reduced moments
    """
    rm = m[:]
    f = 1.0
    for i in range(len(m)):
        f = f / (i+1)
        rm[i] *= f
    return rm

def MomsFromReducedMoms (rm):
    """
    Returns the raw moments given the reduced moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The reduced moments are: `\displaystyle r_i=\frac{m_i}{i!}`
       
    Parameters
    ----------
    rm : vector of doubles
        The list of reduced moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments
    """
    m = rm[:]
    f = 1.0
    for i in range(len(m)):
        f = f * (i+1)
        m[i] *= f
    return m

def MomsFromFactorialMoms (fm):
    """
    Returns the raw moments given the factorial moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The factorial moments are: `f_i=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1))`
       
    Parameters
    ----------
    fm : vector of doubles
        The list of factorial moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments
    
    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    
    """
    n = len(fm)
    m = [fm[0]]
    for i in range (1,n):
        eh = -np.poly(range(i+1))[i:0:-1]
        m.append(fm[i]+np.dot(eh,m[0:i]))
    return m
    
def FactorialMomsFromMoms (m):   
    """
    Returns the factorial moments given the raw moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The factorial moments are: `f_i=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1))`
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    fm : vector of doubles
        The list of factorial moments
    
    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    
    """
    n = len(m)
    fm = []
    for i in range(n):
        eh = np.poly(range(i+1))[i::-1]
        fm.append (np.dot(eh, m[0:i+1]))
    return fm
    
def HankelMomsFromMoms (m):
    """
    Returns the Hankel moments given the raw moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The ith Hankel moment is the determinant of matrix 
    `\Delta_{i/2}`, if i is even, 
    and it is the determinant of `\Delta^{(1)}_{(i+1)/2}`, 
    if i is odd. For the definition of matrices `\Delta`
    and `\Delta^{(1)}` see [1]_.
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    hm : vector of doubles
        The list of Hankel moments
    
    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem
    """
    hm = []
    for i in range(len(m)):
        if i%2 == 0:
            N = i//2 + 1
            H = la.hankel(m[0:N], m[N-1:2*N-1])
        else:
            N = (i+1) // 2 + 1
            H = la.hankel([1] + m[0:N-1], m[N-2:2*N-2])
        hm.append (la.det(H))
    return hm
    
def MomsFromHankelMoms (hm):
    """
    Returns the raw moments given the Hankel moments.
    
    The raw moments are: `m_i=E(\mathcal{X}^i)`
    
    The ith Hankel moment is the determinant of matrix 
    `\Delta_{i/2}`, if i is even, 
    and it is the determinant of `\Delta^{(1)}_{(i+1)/2}`, 
    if i is odd. For the definition of matrices `\Delta`
    and `\Delta^{(1)}` see [1]_.
       
    Parameters
    ----------
    hm : vector of doubles
        The list of Hankel moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments
    
    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem
    """
    m = [hm[0]]
    for i in range(1,len(hm)):
        if i%2 == 0:
            N = i//2 + 1
            H = la.hankel(m[0:N], m[N-1:2*N-2] + [0])
        else:
            N = (i+1) // 2 + 1
            H = la.hankel([1] + m[0:N-1], m[N-2:2*N-3] + [0])
        h = hm[i]
        rH = np.delete(H, N-1, 0)
        for j in range(N):
            cofactor = (-1)**(N+j-1) * la.det(np.delete(rH, j, 1))
            if j<N-1:
                h -= cofactor * H[N-1,j]
            else:
                m.append(h/cofactor)
    return m

def JMomsFromJFactorialMoms (jfmoms):
    """
    Returns the lag-1 joint raw moments given the lag-1 joint
    factorial moments.
    
    The lag-1 joint raw moments are: 
    `m_{i,j}=E(\mathcal{X}^i\mathcal{Y}^j)`
    
    The factorial moments are: 
    `f_{ij}=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1)\mathcal{Y}(\mathcal{Y}-1)\cdots(\mathcal{Y}-j+1))`
       
    Parameters
    ----------
    jfm : matrix, shape (M,M)
        The matrix of joint factorial moments. The entry in 
        row i and column j is `f_{i,j},i\geq 1,j\geq 1`.
        
    Returns
    -------
    jm : matrix, shape (M,M)
        The matrix of joint raw moments. The entry in row i
        and column j is `m_{i,j},i\geq 1,j\geq 1`.
        
    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    
    """
    jfmoms = np.array(jfmoms)
    s1 = jfmoms.shape[0]
    s2 = jfmoms.shape[1]
    jmoms = np.zeros((s1,s2))
    for i in range(s1):
        for j in range(s2):
            xCoeff = np.poly(range(i+1))[i::-1]
            yCoeff = np.poly(range(j+1))[j::-1]
            eh = -np.outer(xCoeff,yCoeff).T
            jmoms[i,j] = jfmoms[i,j] + np.trace(jmoms[:i+1,:j+1].dot(eh))
    return jmoms

def JFactorialMomsFromJMoms (jmoms):
    """
    Returns the lag-1 joint factorial moments given the 
    lag-1 joint raw moments.
    
    The lag-1 joint raw moments are: 
    `m_{i,j}=E(\mathcal{X}^i\mathcal{Y}^j)`
    
    The factorial moments are: 
    `f_{ij}=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1)\mathcal{Y}(\mathcal{Y}-1)\cdots(\mathcal{Y}-j+1))`
       
    Parameters
    ----------
    jm : matrix, shape (M,M)
        The matrix of joint raw moments. The entry in row i
        and column j is `m_{i,j},i\geq 1,j\geq 1`.
        
    Returns
    -------
    jfm : matrix, shape (M,M)
        The matrix of joint factorial moments. The entry in 
        row i and column j is `f_{i,j},i\geq 1,j\geq 1`.
    
    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    
    """
    jmoms = np.array(jmoms)
    s1 = jmoms.shape[0]
    s2 = jmoms.shape[1]
    jfmoms = np.zeros((s1,s2))
    for i in range(s1):
        for j in range(s2):
            xCoeff = np.poly(range(i+1))[i::-1]
            yCoeff = np.poly(range(j+1))[j::-1]
            eh = np.outer(xCoeff,yCoeff).T
            jfmoms[i,j] = np.trace(jmoms[:i+1,:j+1].dot(eh))
    return jfmoms
