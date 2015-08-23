butools.dmap.DMMAPFromDMRAP
===========================

.. currentmodule:: butools.dmap

.. np:function:: DMMAPFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = DMMAPFromDMRAP(H, precision)`
        * - Mathematica:
          - :code:`D = DMMAPFromDMRAP[H, precision]`
        * - Python/Numpy:
          - :code:`D = DMMAPFromDMRAP(H, precision)`

    Obtains a Markovian representation of a discrete rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to transform
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP (if found)

    References
    ----------
    .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
           traffic based decomposition of two-class queueing 
           networks with priority service". COMPUTER NETWORKS 
           53:(8) pp. 1235-1248. (2009)

    Examples
    ========
    For Matlab:

    >>> H0 = [0.15, 0.2, 0.18; -0.20, 0.17, 0.22; 0.19, 0.15, 0.16];
    >>> H1 = [0.01, 0.08, 0.16; 0.02, 0.2, 0.07; 0.02, 0.15, 0.17];
    >>> H2 = [0.14, 0.07, 0.01; 0.19, 0.02, 0.31; 0.06, 0.1, 0.];
    >>> H = {H0, H1, H2};
    >>> moms = MarginalMomentsFromDMRAP(H);
    >>> disp(moms);
           1.6264       3.6055       10.991       43.903       218.08
    >>> jmom = LagkJointMomentsFromDMRAP(H, 3, 1);
    >>> G = DMMAPFromDMRAP(H);
    >>> disp(G{1});
          0.12149      0.28833      0.11968
       5.6441e-06      0.17495     0.068392
          0.14667      0.12596      0.18355
    >>> disp(G{2});
         0.026095    5.788e-06      0.13053
         0.069062      0.26926   1.1436e-05
         0.073725      0.21205     0.084643
    >>> disp(G{3});
          0.14939    0.0019471      0.16253
          0.12377     0.010576      0.28397
         0.047622      0.12573   3.4639e-05
    >>> rmoms = MarginalMomentsFromDMMAP(G);
    >>> disp(rmoms);
           1.6264       3.6055       10.991       43.903       218.08
    >>> rjmom = LagkJointMomentsFromDMMAP(G, 3, 1);
    >>> err = norm(rjmom{1}-jmom{1})+norm(rjmom{2}-jmom{2});
    >>> disp(err);
       7.8933e-13

    For Mathematica:

    
    For Python/Numpy:

    >>> H0 = ml.matrix([[0.15, 0.2, 0.18],[-0.20, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1 = ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2 = ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.31],[0.06, 0.1, 0.]])
    >>> H = [H0, H1, H2]
    >>> moms = MarginalMomentsFromDMRAP(H)
    >>> print(moms)
    [1.6263896740154515, 3.6054695734649633, 10.991320699229343, 43.902870881249427, 218.07910677758866]
    >>> jmom = LagkJointMomentsFromDMRAP(H, 3, 1)
    >>> G = DMMAPFromDMRAP(H)
    >>> print(G[0])
    [[ 0.15737  0.31966  0.07273]
     [ 0.01095  0.22961  0.06546]
     [ 0.16355  0.1105   0.09302]]
    >>> print(G[1])
    [[ 0.08443  0.03798  0.09227]
     [ 0.01529  0.2798   0.01004]
     [ 0.02635  0.39383  0.01577]]
    >>> print(G[2])
    [[ 0.12909  0.09375  0.01271]
     [ 0.22947  0.01525  0.14412]
     [ 0.01002  0.1713   0.01567]]
    >>> rmoms = MarginalMomentsFromDMMAP(G)
    >>> print(rmoms)
    [1.6263896740154513, 3.6054695734649629, 10.991320699229343, 43.902870881249434, 218.07910677758866]
    >>> rjmom = LagkJointMomentsFromDMMAP(G, 3, 1)
    >>> err = la.norm(rjmom[0]-jmom[0])+la.norm(rjmom[1]-jmom[1])
    >>> print(err)
    6.76975985004e-14

