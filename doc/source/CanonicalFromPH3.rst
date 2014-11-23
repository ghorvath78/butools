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
    --------
    For Matlab:

    >>> a=[0.1 0.9 0];
    >>> A=[-6.2 2 0; 2 -9 1; 1 0 -3];
    >>> [b,B]=CanonicalFromPH3(a,A);
    >>> b
      0.58305      0.32736     0.089589
    >>> B
      -9.9819            0            0
       5.3405      -5.3405            0
            0       2.8776      -2.8776
    >>> C=SimilarityMatrix(A,B);    
    >>> norm(A*C-C*B)
      8.4853e-15
    >>> norm(a*C-b)
      8.3463e-16
      
    For Python/Numpy:
    
    >>> a=ml.matrix([[0.1, 0.9, 0]])
    >>> A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> b,B=CanonicalFromPH3(a,A)
    >>> print(b)
    [[ 0.58305413  0.32735661  0.08958926]]
    >>> print(B)
    [[-9.98192063  0.          0.        ]
     [ 5.34047103 -5.34047103  0.        ]
     [ 0.          2.87760834 -2.87760834]]
    >>> C=SimilarityMatrix(A,B)
    >>> print(la.norm(A*C-C*B))
    8.94161968959e-14
    >>> print(la.norm(a*C-b))
    3.02144014162e-15

