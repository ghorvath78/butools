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
    
