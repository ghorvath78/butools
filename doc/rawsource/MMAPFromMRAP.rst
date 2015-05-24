butools.map.MMAPFromMRAP
========================

.. currentmodule:: butools.map

.. np:function:: MMAPFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = MMAPFromMRAP(H, precision)`
        * - Mathematica:
          - :code:`D = MMAPFromMRAP[H, precision]`
        * - Python/Numpy:
          - :code:`D = MMAPFromMRAP(H, precision)`

    Obtains a Markovian representation of a rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP to transform
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP (if found)

    References
    ----------
    .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
           traffic based decomposition of two-class queueing 
           networks with priority service". COMPUTER NETWORKS 
           53:(8) pp. 1235-1248. (2009)

    Examples
    --------
    For Matlab:

    >>> H0=[-5 0.28 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9];
    >>> H1=[-0.08 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1];
    >>> H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8];
    >>> CheckMMAPRepresentation({H0,H1,H2});
    CheckMMAPRepresentation: Some of the matrices H1 ... HM have a negative element!
         0
    >>> D=MMAPFromMRAP(H);
    >>> D{1}
          -4.6311      0.17674      0.85564      0.92969
           1.0899      -8.0182       0.6718     0.037188
           1.2592     0.039062      -4.2279      0.92969
          0.85625       2.0472       3.0031      -9.1228
    >>> D{2}
         0.013162      0.54118      0.10313      0.19567
          0.21755      0.95606       1.8821     0.088018
           0.2338     0.055469     0.059953      0.65078
          0.96289     0.065679     0.030318      0.09082
    >>> D{3}
           0.3392     0.016891       0.1532       1.3066
           2.3529     0.093338      0.62233    0.0070033
          0.14415     0.088281       0.6746     0.092969
          0.41951      0.96644     0.087782      0.59286
    >>> CheckMMAPRepresentation(D)
         1
    >>> jmom=LagkJointMomentsFromMRAP({H0,H1,H2},3,1);
    >>> rjmom=LagkJointMomentsFromMMAP(D,3,1);
    >>> error = norm(rjmom{1}-jmom{1}) + norm(rjmom{2}-jmom{2});
       1.1874e-15

    For Python/Numpy:
    
    >>> H0=ml.matrix([[-5, 0.28, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    >>> H1=ml.matrix([[-0.08, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    >>> H=(H0,H1,H2)
    >>> print(CheckMMAPRepresentation(H))
    CheckMMAPRepresentation: Some of the matrices H1 ... HM have a negative element!
    False
    >>> D=MMAPFromMRAP(H)
    >>> print(D[0])
    [[-4.64977368  0.17674317  0.87431693  0.9296875 ]
     [ 1.0891793  -8.01818988  0.67254576  0.0371875 ]
     [ 1.24047394  0.0390625  -4.20922394  0.9296875 ]
     [ 0.92119055  2.04994681  2.9219339  -9.1228125 ]]
    >>> print(D[1])
    [[ 0.0092309   0.54117557  0.10706401  0.19567014]
     [ 0.21578177  0.95606452  1.88387129  0.08801829]
     [ 0.22072296  0.05546875  0.07302704  0.65078125]
     [ 0.956814    0.07543657  0.0328261   0.08167755]]
    >>> print(D[2])
    [[ 0.31294725  0.01689127  0.17944857  1.30659838]
     [ 2.35272587  0.09333804  0.62247422  0.00700332]
     [ 0.14228058  0.08828125  0.67646942  0.09296875]
     [ 0.41102999  0.96500448  0.08970727  0.61724529]]
    >>> print(CheckMMAPRepresentation(D))
    True
    >>> rmoms=MarginalMomentsFromMMAP(D)
    >>> rjmom=LagkJointMomentsFromMMAP(D,3,1)
    >>> error = la.norm(rjmom[0]-jmom[0]) + la.norm(rjmom[1]-jmom[1])
    >>> print(error)
    9.65812159364e-16

