butools.dph.CanonicalFromDPH2
=============================

.. currentmodule:: butools.dph

.. np:function:: CanonicalFromDPH2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = CanonicalFromDPH2(alpha, A, prec)`
        * - Mathematica:
          - :code:`{beta, B} = CanonicalFromDPH2[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`beta, B = CanonicalFromDPH2(alpha, A, prec)`

    Returns the canonical form of an order-2 discrete phase-type 
    distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,2)
        Initial vector of the discrete phase-type distribution
    A : matrix, shape (2,2)
        Transition probability matrix of the discrete phase-type
        distribution
    prec : double, optional
      Numerical precision for checking the input, default value
      is 1e-14

    Returns
    -------
    beta : matrix, shape (1,2)
      The initial probability vector of the canonical form
    B : matrix, shape (2,2)
      Transition probability matrix of the canonical form

    Examples
    ========
    For Matlab:

    >>> a = [0,1.0];
    >>> A = [0.23, 0.22; 0.41, 0.48];
    >>> [b, B] = CanonicalFromDPH2(a, A);
    >>> disp(b);
          0.88663      0.11337
    >>> disp(B);
          0.68031      0.31969
                0     0.029692
    >>> ev = eig(A);
    >>> disp(ev);
         0.029692
          0.68031
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> a = [1.0,0];
    >>> A = [0, 0.61; 0.56, 0.44];
    >>> [b, B] = CanonicalFromDPH2(a, A);
    >>> disp(b);
      -5.0834e-16            1
    >>> disp(B);
             0.44         0.56
             0.61            0
    >>> ev = eig(A);
    >>> disp(ev);
          -0.4045
           0.8445
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
       4.4871e-16

    For Mathematica:

    >>> a = {0,1.0};
    >>> A = {{0.23, 0.22},{0.41, 0.48}};
    >>> {b, B} = CanonicalFromDPH2[a, A];
    >>> Print[b];
    {0.8866338818412278, 0.1133661181587723}
    >>> Print[B];
    {{0.6803075467922624, 0.31969245320773765},
     {0, 0.029692453207737557}}
    >>> ev = Eigenvalues[A];
    >>> Print[ev];
    {0.6803075467922624, 0.029692453207737557}
    >>> flag = CheckDPHRepresentation[b, B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A, B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> a = {1.0,0};
    >>> A = {{0, 0.61},{0.56, 0.44}};
    >>> {b, B} = CanonicalFromDPH2[a, A];
    >>> Print[b];
    {-5.083438757441196*^-16, 1.0000000000000002}
    >>> Print[B];
    {{0.4400000000000001, 0.5599999999999998},
     {0.6100000000000002, 0}}
    >>> ev = Eigenvalues[A];
    >>> Print[ev];
    {0.8444997998398399, -0.4044997998398398}
    >>> flag = CheckDPHRepresentation[b, B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A, B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1, err2]];
    1.059520934796808*^-15

    For Python/Numpy:

    >>> a = ml.matrix([[0,1.0]])
    >>> A = ml.matrix([[0.23, 0.22],[0.41, 0.48]])
    >>> b, B = CanonicalFromDPH2(a, A)
    >>> print(b)
    [[ 0.88663  0.11337]]
    >>> print(B)
    [[ 0.68031  0.31969]
     [ 0.       0.02969]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [ 0.02969+0.j  0.68031+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> a = ml.matrix([[1.0,0]])
    >>> A = ml.matrix([[0, 0.61],[0.56, 0.44]])
    >>> b, B = CanonicalFromDPH2(a, A)
    >>> print(b)
    [[ -5.08344e-16   1.00000e+00]]
    >>> print(B)
    [[ 0.44  0.56]
     [ 0.61  0.  ]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [-0.4045+0.j  0.8445+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    3.90921887111e-16

