butools.mam.GM1StationaryDistr
==============================

.. currentmodule:: butools.mam

.. np:function:: GM1StationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = GM1StationaryDistr (B, R, K)`
        * - Mathematica:
          - :code:`pi = GM1StationaryDistr [B, R, K]`
        * - Python/Numpy:
          - :code:`pi = GM1StationaryDistr (B, R, K)`

    Returns the stationary distribution of the G/M/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : matrix, shape (N,M*N)
        Matrix blocks of the G/M/1 type generator in the
        regular part, from 0 to M-1, concatenated 
        horizontally.
    B : matrix, shape (N,M*N)
        Matrix blocks of the G/M/1 type generator at the
        boundary, from 0 to M-1, concatenated horizontally.
    R : matrix, shape (N,N)
        Matrix R of the G/M/1 type Markov chain
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
    
    >>> A0 = [0.1, 0; 0, 0.1];
    >>> A1 = [0, 0.2; 0, 0.2];
    >>> A2 = [0, 0.1; 0, 0];
    >>> A3 = [0.3, 0.2; 0.3, 0.2];
    >>> A4 = [0, 0.1; 0.2, 0];
    >>> B0 = [0.7, 0.2; 0.3, 0.6];
    >>> B1 = [0.3, 0.4; 0.5, 0.2];
    >>> B2 = [0.2, 0.4; 0.1, 0.6];
    >>> B3 = [0, 0.1; 0.2, 0];
    >>> A = [A0,A1,A2,A3,A4];
    >>> B = [B0,B1,B2,B3];
    >>> R = GM1FundamentalMatrix(A);
    >>> pi = GM1StationaryDistr(B,R,3);
    >>> pi(1:2)
       0.5152      0.35773
    >>> pi(3:4)
      0.05209     0.058853
    >>> pi(5:6)
        0.0052815    0.0088014

    For Python/Numpy:
       
    >>> B0 = ml.matrix([[0.7, 0.2],[0.3, 0.6]])
    >>> B1 = ml.matrix([[0.3, 0.4],[0.5, 0.2]])
    >>> B2 = ml.matrix([[0.2, 0.4],[0.1, 0.6]])
    >>> B3 = ml.matrix([[0, 0.1],[0.2, 0]])
    >>> A0 = ml.matrix([[0.1, 0],[0, 0.1]])
    >>> A1 = ml.matrix([[0, 0.2],[0, 0.2]])
    >>> A2 = ml.matrix([[0, 0.1],[0, 0]])
    >>> A3 = ml.matrix([[0.3, 0.2],[0.3, 0.2]])
    >>> A4 = ml.matrix([[0, 0.1],[0.2, 0]])
    >>> B = (B0,B1,B2,B3)
    >>> A = (A0,A1,A2,A3,A4)
    >>> R = GM1FundamentalMatrix (A)
    >>> pi = GM1StationaryDistr (B,R,3)
    >>> print(pi[0,0:2])
    [[ 0.51519712  0.35773254]]
    >>> print(pi[0,2:4])
    [[ 0.05208979  0.05885256]]
    >>> print(pi[0,4:6])
    [[ 0.00528148  0.0088014 ]]

