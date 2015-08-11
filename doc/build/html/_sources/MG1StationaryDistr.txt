butools.mam.MG1StationaryDistr
==============================

.. currentmodule:: butools.mam

.. np:function:: MG1StationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = MG1StationaryDistr (A, B, G, K)`
        * - Mathematica:
          - :code:`pi = MG1StationaryDistr [A, B, G, K]`
        * - Python/Numpy:
          - :code:`pi = MG1StationaryDistr (A, B, G, K)`

    Returns the stationary distribution of the M/G/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : matrix, shape (N,M*N)
        Matrix blocks of the M/G/1 type generator in the
        regular part, from 0 to M-1, concatenated 
        horizontally.
    B : matrix, shape (N,M*N)
        Matrix blocks of the M/G/1 type generator at the
        boundary, from 0 to M-1, concatenated horizontally.
    G : matrix, shape (N,N)
        Matrix G of the M/G/1 type Markov chain
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : matrix, shape (1,(K+1)*N)
        The stationary probability vector up to level K
    
    Examples
    --------
    For Matlab:
    
    >>> B0 = [0.1, 0.5; 0.3, 0.4];
    >>> B1 = [0, 0.1; 0, 0];
    >>> B2 = [0.2, 0; 0, 0.2];
    >>> B3 = [0, 0.1; 0.1, 0];
    >>> A0 = [0.4, 0.2; 0.3, 0.4];
    >>> A1 = [0, 0.1; 0, 0];
    >>> A2 = [0, 0.2; 0, 0.2];
    >>> A3 = [0.1, 0; 0.1, 0];
    >>> B = [B0,B1,B2,B3];
    >>> A = [A0,A1,A2,A3];
    >>> G = MG1FundamentalMatrix (A);
    >>> pi = MG1StationaryDistr (A,B,G,3);
    >>> pi(1:2)
          0.10293      0.15492
    >>> pi(3:4)
         0.060187      0.07362
    >>> pi(5:6)
         0.068822      0.10886

    For Python/Numpy:

    >>> B0 = ml.matrix([[0.1, 0.5],[0.3, 0.4]])
    >>> B1 = ml.matrix([[0, 0.1],[0, 0]])
    >>> B2 = ml.matrix([[0.2, 0],[0, 0.2]])
    >>> B3 = ml.matrix([[0, 0.1],[0.1, 0]])
    >>> A0 = ml.matrix([[0.4, 0.2],[0.3, 0.4]])
    >>> A1 = ml.matrix([[0, 0.1],[0, 0]])
    >>> A2 = ml.matrix([[0, 0.2],[0, 0.2]])
    >>> A3 = ml.matrix([[0.1, 0],[0.1, 0]])
    >>> B = (B0,B1,B2,B3)
    >>> A = (A0,A1,A2,A3)
    >>> G = MG1FundamentalMatrix(A)
    >>> pi = MG1StationaryDistr (A,B,G,3)
    >>> print(pi[0,0:2])
    [[ 0.10292778  0.15491506]]
    >>> print(pi[0,2:4])
    [[ 0.06018653  0.07361959]]
    >>> print(pi[0,4:6])
    [[ 0.06882166  0.10885956]]

