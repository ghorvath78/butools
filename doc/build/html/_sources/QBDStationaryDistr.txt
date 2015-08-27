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
    pi : array, length (K+1)*N
        The stationary probability vector up to level K
    
    Examples
    ========
    For Matlab:

    >>> B = [0., 0.; 3., 4.];
    >>> L = [-6., 5.; 3., -12.];
    >>> F = [1., 0.; 2., 0.];
    >>> L0 = [-6., 5.; 6., -8.];
    >>> pi = QBDStationaryDistr(pi0, R, 5);
    >>> disp(pi);
      Columns 1 through 6
          0.22992      0.18681      0.16802     0.086221     0.094781     0.048638
      Columns 7 through 12
         0.053466     0.027437     0.030161     0.015477     0.017014    0.0087307

    For Mathematica:

    >>> B = {{0., 0.},{3., 4.}};
    >>> L = {{-6., 5.},{3., -12.}};
    >>> F = {{1., 0.},{2., 0.}};
    >>> L0 = {{-6., 5.},{6., -8.}};
    >>> pi = QBDStationaryDistr[pi0, R, 5];
    >>> Print[pi];
    {0.22992392223161465, 0.18681318681318687, 0.1680213277846414, 0.08622147083685547, 0.09478126182723359, 0.048637752779764606, 0.05346635282561894, 0.02743668105525182, 0.03016050672214401, 0.015477102133731794, 0.017013619176594053, 0.00873067299851537}

    For Python/Numpy:

    >>> B = ml.matrix([[0., 0.],[3., 4.]])
    >>> L = ml.matrix([[-6., 5.],[3., -12.]])
    >>> F = ml.matrix([[1., 0.],[2., 0.]])
    >>> L0 = ml.matrix([[-6., 5.],[6., -8.]])
    >>> pi = QBDStationaryDistr(pi0, R, 5)
    >>> print(pi)
    [[ 0.22992  0.18681  0.16802  0.08622  0.09478  0.04864  0.05347  0.02744  0.03016  0.01548  0.01701  0.00873]]

