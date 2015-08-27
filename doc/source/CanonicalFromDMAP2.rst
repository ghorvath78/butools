butools.dmap.CanonicalFromDMAP2
===============================

.. currentmodule:: butools.dmap

.. np:function:: CanonicalFromDMAP2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[G0, G1] = CanonicalFromDMAP2(D0, D1, prec)`
        * - Mathematica:
          - :code:`{G0, G1} = CanonicalFromDMAP2[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`G0, G1 = CanonicalFromDMAP2(D0, D1, prec)`

    Returns the canonical form of an order-2 discrete Markovian
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (2,2)
        The D0 matrix of the DMAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the DMAP(2)
    prec : double, optional
        Numerical precision to check the input, default 
        value is 1e-14

    Returns
    -------
    G0 : matrix, shape (1,2)
        The D0 matrix of the canonical DMAP(2)
    G1 : matrix, shape (2,2)
        The D1 matrix of the canonical DMAP(2)

    Examples
    ========
    For Matlab:

    >>> D0 = [0.46, 0.28; 0.35, 0.23];
    >>> D1 = [0.08, 0.18; 0.14, 0.28];
    >>> [H0, H1] = CanonicalFromDMAP2(D0, D1);
    >>> disp(H0);
           0.6785      0.31704
                0     0.011496
    >>> disp(H1);
                0     0.004455
           0.6285         0.36
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
       8.9731e-14
    >>> D0 = [0.26, 0.28; 0.35, 0.23];
    >>> D1 = [0.28, 0.18; 0.14, 0.28];
    >>> [H0, H1] = CanonicalFromDMAP2(D0, D1);
    >>> disp(H0);
             0.49      0.38875
         0.098265            0
    >>> disp(H1);
          0.12125            0
          0.46299      0.43875
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
       1.2787e-14
    >>> D0 = [0.14, 0.34; 0.35, 0.23];
    >>> D1 = [0.22, 0.3; 0.28, 0.14];
    >>> [H0, H1] = CanonicalFromDMAP2(D0, D1);
    >>> disp(H0);
             0.37      0.51734
          0.16778            0
    >>> disp(H1);
                0      0.11266
          0.47222         0.36
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
       1.2567e-15

    For Mathematica:

    >>> D0 = {{0.46, 0.28},{0.35, 0.23}};
    >>> D1 = {{0.08, 0.18},{0.14, 0.28}};
    >>> {H0, H1} = CanonicalFromDMAP2[D0, D1];
    >>> Print[H0];
    {{0.6785041229130464, 0.31704085460114795},
     {0, 0.011495877086969752}}
    >>> Print[H1];
    {{0, 0.004455022485805612},
     {0.6285041229130266, 0.3600000000000036}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    8.953600966260704*^-14
    >>> D0 = {{0.26, 0.28},{0.35, 0.23}};
    >>> D1 = {{0.28, 0.18},{0.14, 0.28}};
    >>> {H0, H1} = CanonicalFromDMAP2[D0, D1];
    >>> Print[H0];
    {{0.4900000000000001, 0.38874507866387564},
     {0.09826490956822952, 0}}
    >>> Print[H1];
    {{0.12125492133612426, 0},
     {0.462990011767894, 0.4387450786638765}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    2.6634419263318992*^-14
    >>> D0 = {{0.14, 0.34},{0.35, 0.23}};
    >>> D1 = {{0.22, 0.3},{0.28, 0.14}};
    >>> {H0, H1} = CanonicalFromDMAP2[D0, D1];
    >>> Print[H0];
    {{0.37, 0.5173403532281082},
     {0.16778122846668353, 0}}
    >>> Print[H1];
    {{0, 0.11265964677189179},
     {0.4722187715333161, 0.3600000000000004}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    5.0785547211252284*^-15

    For Python/Numpy:

    >>> D0 = ml.matrix([[0.46, 0.28],[0.35, 0.23]])
    >>> D1 = ml.matrix([[0.08, 0.18],[0.14, 0.28]])
    >>> H0, H1 = CanonicalFromDMAP2(D0, D1)
    >>> print(H0)
    [[ 0.6785   0.31704]
     [ 0.       0.0115 ]]
    >>> print(H1)
    [[ 0.       0.00446]
     [ 0.6285   0.36   ]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    8.9717085547e-14
    >>> D0 = ml.matrix([[0.26, 0.28],[0.35, 0.23]])
    >>> D1 = ml.matrix([[0.28, 0.18],[0.14, 0.28]])
    >>> H0, H1 = CanonicalFromDMAP2(D0, D1)
    >>> print(H0)
    [[ 0.49     0.38875]
     [ 0.09826  0.     ]]
    >>> print(H1)
    [[ 0.12125  0.     ]
     [ 0.46299  0.43875]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    1.74838096756e-15
    >>> D0 = ml.matrix([[0.14, 0.34],[0.35, 0.23]])
    >>> D1 = ml.matrix([[0.22, 0.3],[0.28, 0.14]])
    >>> H0, H1 = CanonicalFromDMAP2(D0, D1)
    >>> print(H0)
    [[ 0.37     0.51734]
     [ 0.16778  0.     ]]
    >>> print(H1)
    [[ 0.       0.11266]
     [ 0.47222  0.36   ]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    1.92296268638e-16

