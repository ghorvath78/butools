butools.map.MinimalRepFromMRAP
==============================

.. currentmodule:: butools.map

.. np:function:: MinimalRepFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = MinimalRepFromMRAP(H, how, precision)`
        * - Mathematica:
          - :code:`H = MinimalRepFromMRAP[H, how, precision]`
        * - Python/Numpy:
          - :code:`H = MinimalRepFromMRAP(H, how, precision)`

    Returns the minimal representation of a marked rational
    arrival process.

    Parameters
    ----------
    H : list of matrices of shape (M,M)
        The list of H0, H1, ..., HK matrices of the marked
        rational arrival process
    how : {"obs", "cont", "obscont"}, optional        
        Determines how the representation is minimized. 
        "cont" means controllability, "obs" means 
        observability, "obscont" means that the rational arrival
        process is minimized in both respects. Default value 
        is "obscont".
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.
    
    Returns
    -------
    D : list of matrices of shape (M,M)
        The D0, D1, ..., DK matrices of the minimal 
        representation

    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.

    Examples
    --------
    For Matlab:
    
    >>> D0=[-5 1 0; 3 -3 0; 1 1 -5];
    >>> D1=[0 0 0.8; 0 0 0; 0.2 0.2 0.2];
    >>> D2=[0 0 3.2; 0 0 0; 0.8 0.8 0.8];
    >>> H=MinimalRepFromMRAP({D0,D1,D2},'obs');    
    >>> H0=H{1}
          -4.4074       1.6931
          0.84259      -2.5926
    >>> H1=H{2}
          0.40741      0.13545
          0.55741     -0.20741
    >>> H2=H{3}
           1.6296       0.5418
           2.2296     -0.82963
    >>> C = SimilarityMatrix(H0,D0);
    >>> dissimilarity = norm(H0*C-C*D0) + norm(H1*C-C*D1) + norm(H2*C-C*D2)
       5.5841e-15
    >>> H=MinimalRepFromMRAP({D0,D1,D2},'cont');    
    >>> H{1}
        -5     1     0
         3    -3     0
         1     1    -5
    >>> H{2}
                0            0          0.8
                0            0            0
              0.2          0.2          0.2
    >>> H{3}
                0            0          3.2
                0            0            0
              0.8          0.8          0.8
    >>> H=MinimalRepFromMRAP({D0,D1,D2},'obscont');
    >>> H{1}
          -4.4074       1.6931
          0.84259      -2.5926
    >>> H{2}
          0.40741      0.13545
          0.55741     -0.20741
    >>> H{3}
           1.6296       0.5418
           2.2296     -0.82963

