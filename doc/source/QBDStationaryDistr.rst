butools.mam.QBDStationaryDistr
==============================

.. currentmodule:: butools.mam

.. np:function:: QBDStationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = QBDStationaryDistr (pi0, R, K)`
        * - Mathematica:
          - :code:`pi = QBDStationaryDistr [pi0, R, K]`
        * - Python/Numpy:
          - :code:`pi = QBDStationaryDistr (pi0, R, K)`

    Returns the stationary distribution of a QBD up to a
    given level K.
    
    Parameters
    ----------
    pi0 : matrix, shape (1,N)
        The stationary probability vector of level zero
    R : matrix, shape (N,N)
        The matrix parameter of the matrix geometrical
        distribution of the QBD 
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : matrix, shape (K+1,N)
        The stationary probability vector up to level K
    
    Examples
    --------
    For Matlab:
    
    >>> B = [0,0;3,4];
    >>> L = [-6,5;3,-12];
    >>> F = [1,0;2,0];
    >>> L0 = [-6,5;6,-8];
    >>> [pi0, R] = QBDSolve (B, L, F, L0);
    >>> pi = QBDStationaryDistr (pi0, R, 5);
    >>> pi(1:2)
          0.22992      0.18681
    >>> pi(3:4)
          0.16802     0.086221
    >>> pi(5:6)
         0.094781     0.048638
    >>> plot(sum(pi,2));

    For Python/Numpy:
    
    >>> B = ml.matrix([[0,0],[3,4]])
    >>> L = ml.matrix([[-6,5],[3,-12]])
    >>> F = ml.matrix([[1,0],[2,0]])
    >>> L0 = ml.matrix([[-6,5],[6,-8]])
    >>> pi0, R = QBDSolve (B, L, F, L0)
    >>> pi = QBDStationaryDistr (pi0, R, 5)
    >>> print(pi[0,0:2])
    [[ 0.22992392  0.18681319]]    
    >>> print(pi[0,2:4])
    [[ 0.16802133  0.08622147]]
    >>> print(pi[0,4:6])    
    [[ 0.09478126  0.04863775]]
    >>> plt.plot(np.sum(pi,1))
    
