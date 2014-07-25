butools.ph.MEOrder
==================

.. currentmodule:: butools.ph

.. np:function:: MEOrder

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`order = MEOrder(alpha, A, kind, prec)`
        * - Mathematica:
          - :code:`order = MEOrder[alpha, A, kind, prec]`
        * - Python/Numpy:
          - :code:`order = MEOrder(alpha, A, kind, prec)`

    Returns the order of the ME distribution (which is not 
    necessarily equal to the size of the representation).

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    kind : {'obs', 'cont', 'obscont', 'moment'}, optional
        Determines which order is computed. Possibilities: 
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
            controllability order,
        'moment': moment order (which is the default).
    prec : double, optional
        Precision used to detect if the determinant of the 
        Hankel matrix is zero (in case of kind="moment" only),
        or the tolerance for the rank calculation. The
        default value is 1e-10.

    Returns
    -------
    order : int
        The order of ME distribution

    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation
            of rational arrival processes." Madrid Conference on
            Qeueuing theory (MCQT), June 2010.

    Examples
    --------
    For Matlab:
    
    >>> a=[1 1 1 1 1 1]/6;
    >>> A=[-1 0 0 0 0 0; 0.5 -2 1 0 0 0; 1 0 -3 1 0 0; 1 0 1 -4 1 0; 4 0 0 0 -5 0; 5 0 0 0 0 -6];
    >>> co=MEOrder(a,A,'cont')
         2    
    >>> oo=MEOrder(a,A,'obs')
         6
    >>> coo=MEOrder(a,A,'obscont')
         2
    >>> mo=MEOrder(a,A,'moment');
         2
     
    >>> a=[2 1]/3;
    >>> A=[-1 1; 0 -3];
    >>> co=MEOrder(a,A,'cont')
         2
    >>> oo=MEOrder(a,A,'obs')
         1
    >>> coo=MEOrder(a,A,'obscont')
         1
    >>> mo=MEOrder(a,A,'moment');
         1
     
    >>> b = [0.2, 0.3, 0.5];
    >>> B = [-1,0,0;0,-3,1;0,-1,-3];
    >>> [a,A] = MonocyclicPHFromME(b,B);
    >>> length(a)
         9
    >>> size(A)
         9     9    
    >>> co=MEOrder(a,A,'cont')
         9
    >>> oo=MEOrder(a,A,'obs')
         3
    >>> coo=MEOrder(a,A,'obscont')
         3
    >>> mo=MEOrder(a,A,'moment');
         3
     
     

