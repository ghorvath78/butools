butools.ph.MinimalRepFromME
===========================

.. currentmodule:: butools.ph

.. np:function:: MinimalRepFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = MinimalRepFromME(alpha, A, how, precision)`
        * - Mathematica:
          - :code:`{beta, B} = MinimalRepFromME[alpha, A, how, precision]`
        * - Python/Numpy:
          - :code:`beta, B = MinimalRepFromME(alpha, A, how, precision)`

    Returns the minimal representation of the given ME 
    distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    how : {"obs", "cont", "obscont", "moment"}, optional        
        Determines how the representation is minimized. 
        Possibilities:
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
            controllability order,
        'moment': moment order (which is the default).
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.

    Returns
    -------
    beta : vector, shape (1,N)
        The initial vector of the minimal representation
    B : matrix, shape (N,N)
        The matrix parameter of the minimal representation

    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation
            of rational arrival processes." Madrid Conference on
            Qeueuing theory (MCQT), June 2010.

    Examples
    --------
    For Matlab:
    
    >>> b = [0.2, 0.3, 0.5];
    >>> B = [-1,0,0;0,-3,1;0,-1,-3];
    >>> [a,A] = MonocyclicPHFromME(b,B)
    >>> length(a)
         9
    >>> size(A)
         9     9    
    >>> [b,B]=MinimalRepFromME(a,A,'cont');
    >>> size(B)
         9     9    
    >>> [b,B]=MinimalRepFromME(a,A,'obs');
    >>> b
            1   1.3878e-17  -1.1102e-16
    >>> B
          -2.8362     0.036222   1.1102e-16
           -16.61      -3.3369       16.042
           1.1643    -0.051724     -0.82688
    >>> [b,B]=MinimalRepFromME(a,A,'obscont');
    >>> b
            1   1.3878e-17  -1.1102e-16
    >>> B
          -2.8362     0.036222   1.1102e-16
           -16.61      -3.3369       16.042
           1.1643    -0.051724     -0.82688
    >>> [b,B]=MinimalRepFromME(a,A,'moment');
    >>> b
          0.33333      0.33333      0.33333
    >>> B
          -2.1905       1.9222      -3.3698
          -1.0769      -2.3906      0.83162
         -0.51037       0.8033      -2.4189
       

