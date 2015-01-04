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
    --------
    For Matlab:
    
    >>> a=[0 1];
    >>> A=[0.23 0.22; 0.41 0.48];
    >>> [b,B]=CanonicalFromDPH2(a,A);
    >>> b
      0.88663      0.11337
    >>> B
      0.68031      0.31969
            0     0.029692
    >>> C=SimilarityMatrix(A,B);
    >>> norm(A*C-C*B)
    1.2477e-16
    >>> norm(a*C-b)
    1.3878e-16
  
    For Python/Numpy:
    
    >>> a=ml.matrix([[0, 1.0]])
    >>> A=ml.matrix([[0.23, 0.22],[0.41, 0.48]])
    >>> b,B=CanonicalFromDPH2(a,A)
    >>> print(b)
    [[ 0.88663388  0.11336612]]   
    >>> print(B)
    [[ 0.68030755  0.31969245]
     [ 0.          0.02969245]]
    >>> C=SimilarityMatrix(A,B)
    >>> print(la.norm(A*C-C*B))
    1.1443916996305594e-16
    >>> print(la.norm(a*C-b))
    1.1857187100668868e-16
        
