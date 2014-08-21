import numpy as np
import scipy.linalg as la

def NormMomsFromMoms (m):
    return [float(m[i])/m[i-1]/m[0] if i>0 else m[0] for i in range(len(m))]

def MomsFromNormMoms (nm):
    m = nm[:]
    for i in range(1,len(nm)):
        m[i] *= m[0] * m[i-1]
    return m

def ReducedMomsFromMoms (m):
    rm = m[:]
    f = 1.0
    for i in range(len(m)):
        f = f / (i+1)
        rm[i] *= f
    return rm

def MomsFromReducedMoms (rm):
    m = rm[:]
    f = 1.0
    for i in range(len(m)):
        f = f * (i+1)
        m[i] *= f
    return m

def MomsFromFactorialMoms (fm):
    n = len(fm)
    m = [fm[0]]
    for i in range (1,n):
        eh = -np.poly(range(i+1))[i:0:-1]
        m.append(fm[i]+np.dot(eh,m[0:i]))
    return m
    
def FactorialMomsFromMoms (m):   
    n = len(m)
    fm = []
    for i in range(n):
        eh = np.poly(range(i+1))[i::-1]
        fm.append (np.dot(eh, m[0:i+1]))
    return fm
    
def HankelMomsFromMoms (m):
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
