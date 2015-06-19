butools.dph.DPHFromMG
=====================

.. currentmodule:: butools.dph

.. np:function:: DPHFromMG

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = DPHFromMG(alpha, A, precision)`
        * - Mathematica:
          - :code:`{beta, B} = DPHFromMG[alpha, A, precision]`
        * - Python/Numpy:
          - :code:`beta, B = DPHFromMG(alpha, A, precision)`

    Obtains a Markovian representation of a matrix 
    geometric distribution of the same size, if possible.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
        distribution.
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer than the precision

    Returns
    -------
    beta : vector, shape (1,M)
        The initial probability vector of the Markovian 
        representation
    B : matrix, shape (M,M)
        Transition probability matrix of the Markovian 
        representation

    References
    ----------
    .. [1] G Horváth, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)

    Examples
    ========
    For Matlab:

    >>> a = [-0.6, 0.3, 1.3];
    >>> A = [0.1, 0.2, 0; 0.3, 0.1, 0.25; -0.3, 0.2, 0.77];
    >>> flag = CheckMGRepresentation(a,A);
    >>> disp(flag);
         1
    >>> flag = CheckDPHRepresentation(a,A);
    CheckProbMatrix: the matrix has negative element (precision: 1e-11)!
    >>> disp(flag);
         0
    >>> [b,B] = DPHFromMG(a,A);
    >>> disp(b);
             0.05       0.1375       0.8125
    >>> disp(B);
              0.1          0.2            0
            0.425      0.06875      0.15625
            0.141      0.01975      0.80125
    >>> flag = CheckDPHRepresentation(b,B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A,B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1,err2));
       4.0704e-16

    For Mathematica:

    >>> a = {-0.6, 0.3, 1.3};
    >>> A = {{0.1, 0.2, 0},{0.3, 0.1, 0.25},{-0.3, 0.2, 0.77}};
    >>> flag = CheckMGRepresentation[a,A];
    >>> Print[flag];
    True
    >>> flag = CheckDPHRepresentation[a,A];
    "CheckProbVector: The vector has negative element!"
    >>> Print[flag];
    False
    >>> {b,B} = DPHFromMG[a,A];
    >>> Print[b];
    {0.050000000000000044, 0.13749999999999998, 0.8125}
    >>> Print[B];
    {{0.1, 0.2, 0.},
     {0.425, 0.06875, 0.15625},
     {0.14100000000000007, 0.019750000000000018, 0.8012500000000001}}
    >>> flag = CheckDPHRepresentation[b,B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A,B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1,err2]];
    4.783309287441108*^-16

    For Python/Numpy:

    >>> a = ml.matrix([[-0.6, 0.3, 1.3]])
    >>> A = ml.matrix([[0.1, 0.2, 0],[0.3, 0.1, 0.25],[-0.3, 0.2, 0.77]])
    >>> flag = CheckMGRepresentation(a,A)
    >>> print(flag)
    True
    >>> flag = CheckDPHRepresentation(a,A)
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> b,B = DPHFromMG(a,A)
    >>> print(b)
    [[ 0.05    0.1375  0.8125]]
    >>> print(B)
    [[ 0.1      0.2      0.     ]
     [ 0.425    0.06875  0.15625]
     [ 0.141    0.01975  0.80125]]
    >>> flag = CheckDPHRepresentation(b,B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A,B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1,err2))
    3.34426942151e-16

