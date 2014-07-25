butools.ph.APHFrom3Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: APHFrom3Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = APHFrom3Moments(moms, maxSize)`
        * - Mathematica:
          - :code:`{alpha, A} = APHFrom3Moments[moms, maxSize]`
        * - Python/Numpy:
          - :code:`alpha, A = APHFrom3Moments(moms, maxSize)`

    Returns an acyclic PH which has the same 3 moments as
    given. If detects the order and the structure 
    automatically to match the given moments.

    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    maxSize : int, optional
      The maximal size of the resulting APH. The default value
      is 100.

    Returns
    -------
    alpha : vector, shape (1,M)
      The initial probability vector of the APH
    A : matrix, shape (M,M)
      Transient generator matrix of the APH
    
    Raises an error if the moments are not feasible with an
    APH of size "maxSize".
    
    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.

    Examples
    --------
    For Matlab:
    
    >>> moms = [10,125,8400];
    >>> [a,A]=APHFrom3Moments(moms);
    >>> a
       1.3212e-05      0.99999            0            0            0            0   
    >>> A
       -0.0022936    0.0022936            0            0            0            0
                0     -0.50029      0.50029            0            0            0
                0            0     -0.50029      0.50029            0            0
                0            0            0     -0.50029      0.50029            0
                0            0            0            0     -0.50029      0.50029
                0            0            0            0            0     -0.50029
    >>> MomentsFromPH(a,A,3);
           10          125         8400

    >>> moms = [10,525,31400];
    >>> [a,A]=APHFrom3Moments(moms);
    >>> a
      0.21179            0            0            0            0            0            0      0.78821
    >>> A
     -0.15079      0.15079            0            0            0            0            0            0
            0     -0.15079      0.15079            0            0            0            0            0
            0            0     -0.15079      0.15079            0            0            0            0
            0            0            0     -0.15079      0.15079            0            0            0
            0            0            0            0     -0.15079      0.15079            0            0
            0            0            0            0            0     -0.15079      0.15079            0
            0            0            0            0            0            0     -0.15079      0.15079
            0            0            0            0            0            0            0      -5.9502
    >>> MomentsFromPH(a,A,3);
           10          525        31400

