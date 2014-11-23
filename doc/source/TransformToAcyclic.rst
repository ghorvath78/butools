butools.reptrans.TransformToAcyclic
===================================

.. currentmodule:: butools.reptrans

.. np:function:: TransformToAcyclic

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = TransformToAcyclic(A, maxSize, precision)`
        * - Mathematica:
          - :code:`B = TransformToAcyclic[A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`B = TransformToAcyclic(A, maxSize, precision)`

    Transforms an arbitrary matrix to a Markovian bi-diagonal 
    matrix.
    
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
    B : matrix, shape (N,N)
        Transient (bi-diagonal) generator matrix of the
        Markovian acyclic representation.
        
    Notes
    -----
    Calls the 'TransformToMonocyclic' procedure if all the 
    eigenvalues are real, otherwise it raises an error if no
    Markovian acyclic generator has been found up to the 
    given size.

    Raises an error if no Markovian acyclic generator 
    has been found up to the given size.
   
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)

    Examples
    --------
    For Matlab:
    
    >>> A = [-0.8 0.8 0; 0.1 -0.3 0.1; 0.2 0 -0.5];
    >>> B=TransformToAcyclic(A);
     -0.1203       0.1203            0
           0      -0.6158       0.6158
           0            0      -0.8639
    >>> C=SimilarityMatrix(A,B);
    >>> A*C
      -0.10221    0.0096089     0.092604
     -0.086842    -0.013158  -2.0936e-09
     -0.053839    -0.072529     -0.17363
    >>> C*B
      -0.10221    0.0096089     0.092604
     -0.086842    -0.013158  -1.9665e-09
     -0.053839    -0.072529     -0.17363

    For Python/Numpy:
    
    >>> A = ml.matrix([[-0.8, 0.8, 0], [0.1, -0.3, 0.1], [0.2, 0, -0.5]])
    >>> B=TransformToAcyclic(A)
    >>> print(B)
    [[-0.12030368  0.12030368  0.        ]
     [ 0.         -0.61579897  0.61579897]
     [ 0.          0.         -0.86389737]]
    >>> C=SimilarityMatrix(A,B)
    >>> print(A*C)
    [[ -1.02212451e-01   9.60885889e-03   9.26035918e-02]
     [ -8.68417847e-02  -1.31582146e-02  -6.95794991e-10]
     [ -5.38390510e-02  -7.25292122e-02  -1.73631737e-01]]
    >>> print(C*B)
    [[ -1.02212457e-01   9.60885905e-03   9.26035917e-02]
     [ -8.68417900e-02  -1.31582145e-02  -6.50973442e-10]
     [ -5.38390542e-02  -7.25292125e-02  -1.73631737e-01]]
    
    
