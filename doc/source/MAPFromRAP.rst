butools.map.MAPFromRAP
======================

.. currentmodule:: butools.map

.. np:function:: MAPFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = MAPFromRAP(H0, H1, precision)`
        * - Mathematica:
          - :code:`{D0, D1} = MAPFromRAP[H0, H1, precision]`
        * - Python/Numpy:
          - :code:`D0, D1 = MAPFromRAP(H0, H1, precision)`

    Obtains a Markovian representation of a rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision

    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process

    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       

    Examples
    --------
    For Matlab:
    
    >>> H0=[-2.4 2; 2 -9];
    >>> H1=[-1.6 2; 3 4];
    >>> CheckMAPRepresentation(H0,H1)
    CheckMAPRepresentation: D1 has negative element (precision: 1e-14)!
         0
    >>> [D0,D1]=MAPFromRAP(H0,H1);
    >>> D0
          -1.8414     0.079468
         0.012509      -9.5586
    >>> D1
         0.024509       1.7374
           7.1706       2.3755
    >>> CheckMAPRepresentation(D0,D1)
         1
    >>> error = norm(LagkJointMomentsFromRAP(H0,H1,3,1)-LagkJointMomentsFromMAP(D0,D1,3,1))
       6.4694e-16
    
