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
    --------
    For Matlab:

    >>> a=[-0.6 0.3 1.3];
    >>> A=[0.1 0.2 0; 0.3 0.1 0.25; -0.3 0.2 0.77];
    >>> CheckDPHRepresentation(a,A)
         0
    >>> [b,B]=DPHFromMG(a,A);
    >>> CheckDPHRepresentation(b,B)
         1
    >>> b
             0.05       0.1375       0.8125
    >>> B
              0.1          0.2            0
            0.425      0.06875      0.15625
            0.141      0.01975      0.80125
    >>> C=SimilarityMatrix(A,B);    
    >>> norm(A*C-C*B)
       2.0713e-16
    >>> norm(a*C-b)
       6.8807e-16
      
