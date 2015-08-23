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
    ========
    For Matlab:

    >>> A = [-0.8, 0.8, 0.; 0.1, -0.3, 0.1; 0.2, 0., -0.5];
    >>> B = TransformToAcyclic(A);
    >>> disp(B);
          -0.1203       0.1203            0
                0      -0.6158       0.6158
                0            0      -0.8639
    >>> Cm = SimilarityMatrix(A, B);
    >>> err = norm(A*Cm-Cm*B);
    >>> disp(err);
       7.0942e-09

    For Mathematica:

    
    For Python/Numpy:

    >>> A = ml.matrix([[-0.8, 0.8, 0.],[0.1, -0.3, 0.1],[0.2, 0., -0.5]])
    >>> B = TransformToAcyclic(A)
    >>> print(B)
    [[-0.1203  0.1203  0.    ]
     [ 0.     -0.6158  0.6158]
     [ 0.      0.     -0.8639]]
    >>> Cm = SimilarityMatrix(A, B)
    >>> err = la.norm(A*Cm-Cm*B)
    >>> print(err)
    8.75449540243e-09

