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
    ========
    For Matlab:

    >>> D0 = [-2., 2.; 2., -9.];
    >>> D1 = [-2., 2.; 3., 4.];
    >>> [H0, H1] = MAPFromRAP(D0, D1);
    >>> disp(H0);
          -1.4699     -0.18514
         -0.04328      -9.5301
    >>> disp(H1);
         -0.18514       1.8401
           7.3883       2.1851
    >>> err = norm(LagkJointMomentsFromRAP(D0, D1, 3, 1)-LagkJointMomentsFromRAP(H0, H1, 3, 1));
    >>> disp(err);
       2.7443e-15
    >>> D0 = [-2.4, 2.; 2., -9.];
    >>> D1 = [-1.6, 2.; 3., 4.];
    >>> [H0, H1] = MAPFromRAP(D0, D1);
    >>> disp(H0);
          -1.8414     0.079468
         0.012509      -9.5586
    >>> disp(H1);
         0.024509       1.7374
           7.1706       2.3755
    >>> err = norm(LagkJointMomentsFromRAP(D0, D1, 3, 1)-LagkJointMomentsFromRAP(H0, H1, 3, 1));
    >>> disp(err);
       6.4694e-16

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-2., 2.],[2., -9.]])
    >>> D1 = ml.matrix([[-2., 2.],[3., 4.]])
    >>> H0, H1 = MAPFromRAP(D0, D1)
    >>> print(H0)
    [[ -1.46887e+00  -3.05193e-02]
     [  4.20579e-15  -9.53113e+00]]
    >>> print(H1)
    [[-0.33715  1.83654]
     [ 7.19398  2.33715]]
    >>> err = la.norm(LagkJointMomentsFromRAP(D0, D1, 3, 1)-LagkJointMomentsFromRAP(H0, H1, 3, 1))
    >>> print(err)
    1.62866154301e-14
    >>> D0 = ml.matrix([[-2.4, 2.],[2., -9.]])
    >>> D1 = ml.matrix([[-1.6, 2.],[3., 4.]])
    >>> H0, H1 = MAPFromRAP(D0, D1)
    >>> print(H0)
    [[-1.84137  0.07947]
     [ 0.01251 -9.55863]]
    >>> print(H1)
    [[ 0.02451  1.7374 ]
     [ 7.17063  2.37549]]
    >>> err = la.norm(LagkJointMomentsFromRAP(D0, D1, 3, 1)-LagkJointMomentsFromRAP(H0, H1, 3, 1))
    >>> print(err)
    4.674749869e-16

