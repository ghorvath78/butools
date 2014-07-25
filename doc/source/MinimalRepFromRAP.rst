butools.map.MinimalRepFromRAP
=============================

.. currentmodule:: butools.map

.. np:function:: MinimalRepFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = MinimalRepFromRAP(H0, H1, how, precision)`
        * - Mathematica:
          - :code:`{D0, D1} = MinimalRepFromRAP[H0, H1, how, precision]`
        * - Python/Numpy:
          - :code:`D0, D1 = MinimalRepFromRAP(H0, H1, how, precision)`

    Returns the minimal representation of a rational arrival 
    process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    how : {"obs", "cont", "obscont"}, optional      
        Determines how the representation is minimized. "cont" 
        means controllability, "obs" means observability, 
        "obscont" means that the rational arrival process is 
        minimized in both respects. The default value is 
        "obscont"
    precision : double, optional
       Precision used by the Staircase algorithm. The default 
       value is 1e-12.
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the minimal representation
    D1 : matrix, shape (M,M)
        The D1 matrix of the minimal representation

    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.

    Examples
    --------
    For Matlab:
    
    >>> D0=[-5 1 0; 3 -3 0; 1 1 -5];
    >>> D1=[0 0 4; 0 0 0; 1 1 1];
    >>> [H0,H1]=MinimalRepFromRAP(D0,D1,'obs');    
    >>> H0
          -4.4074       1.6931
          0.84259      -2.5926
    >>> H1
            2.037      0.67725
            2.787       -1.037
    >>> C = SimilarityMatrix(H0,D0);
    >>> dissimilarity = norm(H0*C-C*D0) + norm(H1*C-C*D1)
       3.0836e-15    
    >>> [H0,H1]=MinimalRepFromRAP(D0,D1,'cont');    
    >>> H0
        -5     1     0
         3    -3     0
         1     1    -5
    >>> H1
         0     0     4
         0     0     0
         1     1     1
    >>> [H0,H1]=MinimalRepFromRAP(D0,D1,'obscont');
    >>> H0
          -4.4074       1.6931
          0.84259      -2.5926
    >>> H1
            2.037      0.67725
            2.787       -1.037
    
    
