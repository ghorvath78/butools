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
    --------
    For Matlab:

    >>> a=[0.12 0.88];
    >>> A=[-1.28 0; 3.94 -3.94];
    >>> [b,B]=CanonicalFromPH2(a,A);
    >>> b
      0.96102     0.038985
    >>> B
        -1.28         1.28
            0        -3.94   
    >>> C=SimilarityMatrix(A,B);    
    >>> norm(A*C-C*B)
      3.8374e-15    
    >>> norm(a*C-b)
      1.8501e-15 
