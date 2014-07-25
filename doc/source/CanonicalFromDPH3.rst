butools.dph.CanonicalFromDPH3
=============================

.. currentmodule:: butools.dph

.. np:function:: CanonicalFromDPH3

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = CanonicalFromDPH3(alpha, A, prec)`
        * - Mathematica:
          - :code:`{beta, B} = CanonicalFromDPH3[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`beta, B = CanonicalFromDPH3(alpha, A, prec)`

    Returns the canonical form of an order-3 discrete phase-type 
    distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,3)
        Initial vector of the discrete phase-type distribution
    A : matrix, shape (3,3)
        Transition probability matrix of the discrete phase-type
        distribution
    prec : double, optional
      Numerical precision for checking the input, default value
      is 1e-14

    Returns
    -------
    beta : matrix, shape (1,3)
      The initial probability vector of the canonical form
    B : matrix, shape (3,3)
      Transition probability matrix of the canonical form

    Examples
    --------
    For Matlab:
    
    >>> a=[0.67 0.07 0.26];
    >>> A=[0.31 0 0; 0.98 0 0.02; 0.88 0.04 0.08];
    >>> [b,B]=CanonicalFromDPH3(a,A);
    >>> b
      0.15814      0.37915       0.4627
    >>> B
         0.31         0.69            0
            0         0.08         0.92
            0   0.00086957            0
    >>> C=SimilarityMatrix(A,B);
    >>> norm(A*C-C*B)
    2.507e-16
    >>> norm(a*C-b)
    7.8505e-17

