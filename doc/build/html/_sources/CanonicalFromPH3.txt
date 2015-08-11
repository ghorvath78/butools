butools.ph.CanonicalFromPH3
===========================

.. currentmodule:: butools.ph

.. np:function:: CanonicalFromPH3

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = CanonicalFromPH3(alpha, A, prec)`
        * - Mathematica:
          - :code:`{beta, B} = CanonicalFromPH3[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`beta, B = CanonicalFromPH3(alpha, A, prec)`

    Returns the canonical form of an order-3 phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,3)
        Initial vector of the phase-type distribution
    A : matrix, shape (3,3)
        Transient generator of the phase-type distribution
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,3)
      The initial probability vector of the canonical form
    B : matrix, shape (3,3)
      Transient generator matrix of the canonical form

    Notes
    -----
    This procedure calculates 5 moments of the input and
    calls 'PH3From5Moments'.

    Examples
    ========
    For Matlab:

    >>> a = [0.1, 0.9, 0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> [b,B] = CanonicalFromPH3(a,A);
    >>> disp(b);
          0.58305      0.32736     0.089589
    >>> disp(B);
          -9.9819            0            0
           5.3405      -5.3405            0
                0       2.8776      -2.8776
    >>> Cm = SimilarityMatrix(A,B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1,err2));
       1.0269e-14

    For Mathematica:

    >>> a = {0.1, 0.9, 0};
    >>> A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}};
    >>> {b,B} = CanonicalFromPH3[a,A];
    >>> Print[b];
    {0.5830541253440302, 0.3273566132692404, 0.08958926138672949}
    >>> Print[B];
    {{-9.981920626264277, 0., 0.},
     {5.340471031780809, -5.340471031780809, 0.},
     {0., 2.8776083419564285, -2.8776083419564285}}
    >>> Cm = SimilarityMatrix[A,B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1,err2]];
    2.7041031044230435*^-13

    For Python/Numpy:

    >>> a = ml.matrix([[0.1, 0.9, 0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> b,B = CanonicalFromPH3(a,A)
    >>> print(b)
    [[ 0.58305  0.32736  0.08959]]
    >>> print(B)
    [[-9.98192  0.       0.     ]
     [ 5.34047 -5.34047  0.     ]
     [ 0.       2.87761 -2.87761]]
    >>> Cm = SimilarityMatrix(A,B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1,err2))
    8.94161968959e-14

