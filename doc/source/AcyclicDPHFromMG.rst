butools.dph.AcyclicDPHFromMG
============================

.. currentmodule:: butools.dph

.. np:function:: AcyclicDPHFromMG

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = AcyclicDPHFromMG(alpha, A, precision)`
        * - Mathematica:
          - :code:`{beta, B} = AcyclicDPHFromMG[alpha, A, precision]`
        * - Python/Numpy:
          - :code:`beta, B = AcyclicDPHFromMG(alpha, A, precision)`

    Transforms a matrix-geometric representation to an acyclic
    DPH representation of the same size, if possible.
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the acyclic discrete
        phase-type representation
    B : matrix, shape (M,M)
        Transition probability matrix of the acyclic discrete
        phase-type representation
    
    Notes
    -----
    Contrary to 'AcyclicPHFromME' of the 'ph' package, this 
    procedure is not able to extend the size in order to obtain
    a Markovian initial vector.

    Raises an error if A has complex eigenvalues. In this case
    the transformation to an acyclic representation is not 
    possible

    Examples
    --------
    For Matlab:
    
    >>> a=[0 0 1]
    >>> A=[0.22 0 0; 0.3 0.1 0.55; 0.26 0 0.73]
    >>> [b,B]=AcyclicDPHFromMG(a,A);
    >>> b
          0.69103      0.29786     0.011111
    >>> B
         0.73         0.27            0
            0         0.22         0.78
            0            0          0.1
    >>> C=SimilarityMatrix(A,B);
    >>> norm(A*C-C*B)
    2.216e-16
    >>> norm(a*C-b)
     0
  

