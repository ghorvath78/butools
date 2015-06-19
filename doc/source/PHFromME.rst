butools.ph.PHFromME
===================

.. currentmodule:: butools.ph

.. np:function:: PHFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = PHFromME(alpha, A, precision)`
        * - Mathematica:
          - :code:`{beta, B} = PHFromME[alpha, A, precision]`
        * - Python/Numpy:
          - :code:`beta, B = PHFromME(alpha, A, precision)`

    Obtains a Markovian representation of a matrix 
    exponential distribution of the same size, if possible.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer than the precision. The default value
        is 1e-14.

    Returns
    -------
    beta : vector, shape (1,M)
        The initial probability vector of the Markovian 
        monocyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        monocyclic representation

    References
    ----------
    .. [1] G Horváth, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)

    Examples
    ========
    For Matlab:

    >>> a = [-0.4, 1.4];
    >>> A = [-3.8, 2; 2, -9];
    >>> flag = CheckMERepresentation(a,A);
    >>> disp(flag);
         1
    >>> flag = CheckPHRepresentation(a,A);
    CheckProbVector: The vector has negative element (precision: 1e-11)!
    >>> disp(flag);
         0
    >>> [b,B] = PHFromME(a,A);
    >>> disp(b);
         0.013037      0.98696
    >>> disp(B);
          -3.2605       2.5924
          0.34843      -9.5395
    >>> flag = CheckPHRepresentation(b,B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A,B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1,err2));
       1.0162e-15
    >>> a = [-0.5, 1.5];
    >>> A = [-3.8, 2; 2, -9];
    >>> flag = CheckMERepresentation(a,A);
    >>> disp(flag);
         1
    >>> flag = CheckPHRepresentation(a,A);
    CheckProbVector: The vector has negative element (precision: 1e-11)!
    >>> disp(flag);
         0
    >>> [b,B] = PHFromME(a,A);
    >>> disp(b);
        0.0057038       0.9943
    >>> disp(B);
          -3.1279       3.0636
         0.017405      -9.6721
    >>> flag = CheckPHRepresentation(b,B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A,B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1,err2));
       3.0445e-15

    For Mathematica:

    >>> a = {-0.4, 1.4};
    >>> A = {{-3.8, 2},{2, -9}};
    >>> flag = CheckMERepresentation[a,A];
    >>> Print[flag];
    True
    >>> flag = CheckPHRepresentation[a,A];
    "CheckProbVector: The vector has negative element!"
    >>> Print[flag];
    False
    >>> {b,B} = PHFromME[a,A];
    >>> Print[b];
    {0.013037109374999953, 0.9869628906249999}
    >>> Print[B];
    {{-3.2604571906887756, 2.5924299798044217},
     {0.3484263627325931, -9.539542809311223}}
    >>> flag = CheckPHRepresentation[b,B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A,B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1,err2]];
    1.1393205654608455*^-15
    >>> a = {-0.5, 1.5};
    >>> A = {{-3.8, 2},{2, -9}};
    >>> flag = CheckMERepresentation[a,A];
    >>> Print[flag];
    True
    >>> flag = CheckPHRepresentation[a,A];
    "CheckProbVector: The vector has negative element!"
    >>> Print[flag];
    False
    >>> {b,B} = PHFromME[a,A];
    >>> Print[b];
    {0.0057037812657654285, 0.9942962187342346}
    >>> Print[B];
    {{-3.1278937744575632, 3.0635853844348873},
     {0.017404720964309853, -9.672106225542432}}
    >>> flag = CheckPHRepresentation[b,B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A,B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1,err2]];
    2.30940485129664*^-15

    For Python/Numpy:

    >>> a = ml.matrix([[-0.4, 1.4]])
    >>> A = ml.matrix([[-3.8, 2],[2, -9]])
    >>> flag = CheckMERepresentation(a,A)
    >>> print(flag)
    True
    >>> flag = CheckPHRepresentation(a,A)
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> b,B = PHFromME(a,A)
    >>> print(b)
    [[ 0.01304  0.98696]]
    >>> print(B)
    [[-3.26046  2.59243]
     [ 0.34843 -9.53954]]
    >>> flag = CheckPHRepresentation(b,B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A,B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1,err2))
    1.18018326364e-15
    >>> a = ml.matrix([[-0.5, 1.5]])
    >>> A = ml.matrix([[-3.8, 2],[2, -9]])
    >>> flag = CheckMERepresentation(a,A)
    >>> print(flag)
    True
    >>> flag = CheckPHRepresentation(a,A)
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> b,B = PHFromME(a,A)
    >>> print(b)
    [[ 0.0057  0.9943]]
    >>> print(B)
    [[-3.12789  3.06359]
     [ 0.0174  -9.67211]]
    >>> flag = CheckPHRepresentation(b,B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A,B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1,err2))
    2.23152066184e-15

