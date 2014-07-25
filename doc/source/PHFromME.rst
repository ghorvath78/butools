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
    --------
    For Matlab:

    >>> a=[-0.4 1.4];
    >>> A=[-3.8 2; 2 -9];
    >>> [b,B]=PHFromME(a,A);
    >>> b
     0.013037      0.98696
    >>> B
      -3.2605       2.5924
      0.34843      -9.5395
    >>> C=SimilarityMatrix(A,B);    
    >>> norm(A*C-C*B)
      1.2949e-15
    >>> norm(a*C-b)
      6.1455e-16
      
