butools.ph.CanonicalFromPH2
===========================

.. currentmodule:: butools.ph

.. np:function:: CanonicalFromPH2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = CanonicalFromPH2(alpha, A, prec)`
        * - Mathematica:
          - :code:`{beta, B} = CanonicalFromPH2[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`beta, B = CanonicalFromPH2(alpha, A, prec)`

    Returns the canonical form of an order-2 phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,2)
      Initial vector of the phase-type distribution
    A : matrix, shape (2,2)
      Transient generator of the phase-type distribution
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    beta : matrix, shape (1,2)
      The initial probability vector of the canonical form
    B : matrix, shape (2,2)
      Transient generator matrix of the canonical form

    Notes
    -----
    This procedure calculates 3 moments of the input and
    calls 'PH2From3Moments'.

    Examples
    ========
    For Matlab:

    >>> a = [0.12,0.88];
    >>> A = [-1.28, 0; 3.94, -3.94];
    >>> [b, B] = CanonicalFromPH2(a, A);
    >>> disp(b);
          0.96102     0.038985
    >>> disp(B);
            -1.28         1.28
                0        -3.94
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
        6.669e-15

    For Mathematica:

    >>> a = {0.12,0.88};
    >>> A = {{-1.28, 0},{3.94, -3.94}};
    >>> {b, B} = CanonicalFromPH2[a, A];
    >>> Print[b];
    {0.9610152284263966, 0.03898477157360336}
    >>> Print[B];
    {{-1.2800000000000014, 1.2800000000000014},
     {0, -3.9399999999999946}}
    >>> Cm = SimilarityMatrix[A, B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1, err2]];
    1.881192080999035*^-15

    For Python/Numpy:

    >>> a = ml.matrix([[0.12,0.88]])
    >>> A = ml.matrix([[-1.28, 0],[3.94, -3.94]])
    >>> b, B = CanonicalFromPH2(a, A)
    >>> print(b)
    [[ 0.96102  0.03898]]
    >>> print(B)
    [[-1.28  1.28]
     [ 0.   -3.94]]
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    4.86095148111e-15

