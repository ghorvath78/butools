butools.moments.CheckMoments
============================

.. currentmodule:: butools.moments

.. np:function:: CheckMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMoments(m, prec)`
        * - Mathematica:
          - :code:`r = CheckMoments[m, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMoments(m, prec)`
    
    Checks if the given moment sequence is valid in the sense
    that it belongs to a distribution with support (0,inf).

    This procedure checks the determinant of :math:`\Delta_n`
    and :math:`\Delta_n^{(1)}` according to [1]_.
    
    Parameters
    ----------
    m : list of doubles, length 2N+1
        The (raw) moments to check 
        (starts with the first moment).
        Its length must be odd.
    prec : double, optional
        Entries with absolute value less than prec are 
        considered to be zeros. The default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem

    Examples
    ========
    For Matlab:

    >>> M = [1.2, 5, 8, 29, 3412];
    >>> flag = CheckMoments(M);
    >>> disp(flag);
         0
    >>> M = [1.3, 2.4, 6.03, 20.5, 89.5];
    >>> flag = CheckMoments(M);
    >>> disp(flag);
         1

    For Mathematica:

    >>> M = {1.2, 5, 8, 29, 3412};
    >>> flag = CheckMoments[M];
    >>> Print[flag];
    False
    >>> M = {1.3, 2.4, 6.03, 20.5, 89.5};
    >>> flag = CheckMoments[M];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> M = [1.2, 5, 8, 29, 3412]
    >>> flag = CheckMoments(M)
    >>> print(flag)
    False
    >>> M = [1.3, 2.4, 6.03, 20.5, 89.5]
    >>> flag = CheckMoments(M)
    >>> print(flag)
    True

