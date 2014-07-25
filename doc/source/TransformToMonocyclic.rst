butools.reptrans.TransformToMonocyclic
======================================

.. currentmodule:: butools.reptrans

.. np:function:: TransformToMonocyclic

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = TransformToMonocyclic(A, maxSize, precision)`
        * - Mathematica:
          - :code:`B = TransformToMonocyclic[A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`B = TransformToMonocyclic(A, maxSize, precision)`

    Transforms an arbitrary matrix to a Markovian monocyclic 
    matrix (see [1]_).
    
    Parameters
    ----------
    A : matrix, shape (N,N)
        Matrix parameter of the initial representation
    maxSize : int, optional
        The maximal order of the resulting Markovian 
        representation. The default value is 100
    precision : double, optional
        Matrix entries smaller than the precision are 
        considered to be zeros. The default value is 1e-14
        
    Returns
    -------
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian monocyclic
        representation. Note that M>N if there are complex 
        eigenvalues.

    Notes
    -----    
    Raises an exception if no Markovian monocyclic generator 
    has been found up to the given size.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)

    Examples
    --------
    For Matlab:
    
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> B=TransformToMonocyclic(A);
       -1        1        0        0        0
        0       -3        3        0        0
        0        0       -3        3        0
        0        0        0       -3        3
        0  0.59259        0        0       -3
    >>> C=SimilarityMatrix(A,B);
    >>> A*C
       1.9399e-15     -0.12308     -0.18462     -0.27692     -0.41538
       6.2774e-17       2.0923      -1.2923      -4.7077       2.9077
       6.6361e-16      0.86154       3.1385      -1.9385      -7.0615
    >>> C*B
       1.9399e-15     -0.12308     -0.18462     -0.27692     -0.41538
       1.1658e-16       2.0923      -1.2923      -4.7077       2.9077
       1.4348e-16      0.86154       3.1385      -1.9385      -7.0615


