butools.map.CanonicalFromMAP2
=============================

.. currentmodule:: butools.map

.. np:function:: CanonicalFromMAP2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[G0, G1] = CanonicalFromMAP2(D0, D1, prec)`
        * - Mathematica:
          - :code:`{G0, G1} = CanonicalFromMAP2[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`G0, G1 = CanonicalFromMAP2(D0, D1, prec)`

    Returns the canonical form of an order-2 Markovian
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (2,2)
        The D0 matrix of the MAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the MAP(2)
    prec : double, optional
        Numerical precision to check the input, default 
        value is 1e-14

    Returns
    -------
    G0 : matrix, shape (1,2)
        The D0 matrix of the canonical MAP(2)
    G1 : matrix, shape (2,2)
        The D1 matrix of the canonical MAP(2)

    Notes
    -----
    This procedure calculates 3 marginal moments and the lag-1
    autocorrelation of the input and calls 'MAP2FromMoments'.

    Examples
    ========
    For Matlab:

    >>> D0 = [-14., 1.; 1., -25.];
    >>> D1 = [6., 7.; 3., 21.];
    >>> [H0, H1] = CanonicalFromMAP2(D0, D1);
    >>> disp(H0);
           -13.91        9.199
                0       -25.09
    >>> disp(H1);
           4.7108            0
            2.801       22.289
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
       2.4485e-13

    For Mathematica:

    >>> D0 = {{-14., 1.},{1., -25.}};
    >>> D1 = {{6., 7.},{3., 21.}};
    >>> {H0, H1} = CanonicalFromMAP2[D0, D1];
    >>> Print[H0];
    {{-13.909830056250456, 9.199027971874015},
     {0, -25.090169943749302}}
    >>> Print[H1];
    {{4.710802084376442, 0},
     {2.8009720281259014, 22.2891979156234}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    2.39754619495158*^-13

    For Python/Numpy:

    >>> D0 = ml.matrix([[-14., 1.],[1., -25.]])
    >>> D1 = ml.matrix([[6., 7.],[3., 21.]])
    >>> H0, H1 = CanonicalFromMAP2(D0, D1)
    >>> print(H0)
    [[-13.90983   9.19903]
     [  0.      -25.09017]]
    >>> print(H1)
    [[  4.7108    0.     ]
     [  2.80097  22.2892 ]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    2.21090873503e-13

