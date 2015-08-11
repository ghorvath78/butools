butools.ph.APH3rdMomentLowerBound
=================================

.. currentmodule:: butools.ph

.. np:function:: APH3rdMomentLowerBound

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m3 = APH3rdMomentLowerBound(m1, m2, n)`
        * - Mathematica:
          - :code:`m3 = APH3rdMomentLowerBound[m1, m2, n]`
        * - Python/Numpy:
          - :code:`m3 = APH3rdMomentLowerBound(m1, m2, n)`

    Returns the lower bound of the third moment of acyclic
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
        The lowest third moment an order-n APH can have with
        the given first and second moment.

    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.

    Examples
    ========
    For Matlab:

    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 3;
    >>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n);
    >>> disp(mom3lower);
           16.577
    >>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n);
    >>> disp(mom3upper);
           17.081
    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 4;
    >>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n);
    >>> disp(mom3lower);
           16.079
    >>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n);
    >>> disp(mom3upper);
       Inf

    For Mathematica:

    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 3;
    >>> mom3lower = APH3rdMomentLowerBound[mean,mom2,n];
    >>> Print[mom3lower];
    16.577458090899114
    >>> mom3upper = APH3rdMomentUpperBound[mean,mom2,n];
    >>> Print[mom3upper];
    17.081405365964272
    >>> mean = 1.9;
    >>> mom2 = 5;
    >>> n = 4;
    >>> mom3lower = APH3rdMomentLowerBound[mean,mom2,n];
    >>> Print[mom3lower];
    16.079377140256383
    >>> mom3upper = APH3rdMomentUpperBound[mean,mom2,n];
    >>> Print[mom3upper];
    Infinity

    For Python/Numpy:

    >>> mean = 1.9
    >>> mom2 = 5
    >>> n = 3
    >>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
    >>> print(mom3lower)
    16.577458090899107
    >>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
    >>> print(mom3upper)
    17.081405365964276
    >>> mean = 1.9
    >>> mom2 = 5
    >>> n = 4
    >>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
    >>> print(mom3lower)
    16.079377140256387
    >>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
    >>> print(mom3upper)
    inf

