butools.ph.APH3rdMomentUpperBound
=================================

.. currentmodule:: butools.ph

.. np:function:: APH3rdMomentUpperBound

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m3 = APH3rdMomentUpperBound(m1, m2, n)`
        * - Mathematica:
          - :code:`m3 = APH3rdMomentUpperBound[m1, m2, n]`
        * - Python/Numpy:
          - :code:`m3 = APH3rdMomentUpperBound(m1, m2, n)`

    Returns the upper bound of the third moment of acyclic
    phase-type (APH) distributions of order n.

    Parameters
    ----------
    m1 : double
        The first moment
    m2 : double
        The second moment
    n : int
        Number of states

    Returns
    -------
    m3 : double
        The highest third moment an order-n APH can have with
        the given first and second moment.

    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.

    Examples
    --------
    For Matlab:
    
    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 3;
    >>> APH3rdMomentLowerBound(mean,mom2,n)
       16.577
    >>> APH3rdMomentUpperBound(mean,mom2,n)
       17.081

    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 4;
    >>> APH3rdMomentLowerBound(mean,mom2,n)
       16.079
    >>> APH3rdMomentUpperBound(mean,mom2,n)
        Inf
 
    For Mathematica:
    
    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 3;
    >>> APH3rdMomentLowerBound[mean,mom2,n]
    16.5775
    >>> APH3rdMomentUpperBound[mean,mom2,n]
    17.0814

    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 4;
    >>> APH3rdMomentLowerBound[mean,mom2,n]
    16.0794
    >>> APH3rdMomentUpperBound[mean,mom2,n]
    Infinity    
 
    For Python/Numpy:
    
    >>> mean = 1.9
    >>> mom2 = 5
    >>> n = 3
    >>> print(APH3rdMomentLowerBound(mean,mom2,n))
    16.577458090899107
    >>> print(APH3rdMomentUpperBound(mean,mom2,n))
    17.081405365964276

    >>> mean = 1.9
    >>> mom2 = 5
    >>> n = 4
    >>> print(APH3rdMomentLowerBound(mean,mom2,n))
    16.079377140256387
    >>> print(APH3rdMomentUpperBound(mean,mom2,n))
    inf
 
