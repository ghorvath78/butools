butools.ph.APH2ndMomentLowerBound
=================================

.. currentmodule:: butools.ph

.. np:function:: APH2ndMomentLowerBound

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m2 = APH2ndMomentLowerBound(m1, n)`
        * - Mathematica:
          - :code:`m2 = APH2ndMomentLowerBound[m1, n]`
        * - Python/Numpy:
          - :code:`m2 = APH2ndMomentLowerBound(m1, n)`

    Returns the lower bound of the second moment of acyclic 
    phase-type (APH) distributions of order n.

    Parameters
    ----------
    m1 : double
        The first moment
    n : int
        Number of states

    Returns
    -------
    m2 : double
        The lowest second moment an order-n APH can have with
        the given first moment.

    References
    ----------
    .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
            discrete and continuous phase-type distributions of
            second order," in In Proc. of UK Performance 
            Evaluation Workshop, UKPEW, 2002"

    Examples
    --------
    For Matlab:
    
    >>> mean = 1.9;
    >>> n = 4;
    >>> mom2 = APH2ndMomentLowerBound(mean,n)
       4.5125
    >>> cv2 = mom2/mean^2-1
         0.25


