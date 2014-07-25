butools.ph.PH2From3Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: PH2From3Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = PH2From3Moments(moms, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = PH2From3Moments[moms, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = PH2From3Moments(moms, prec)`

    Returns a PH(2) which has the same 3 moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,2)
      The initial probability vector of the PH(2)
    A : matrix, shape (2,2)
      Transient generator matrix of the PH(2)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a PH(2).

    References
    ----------
    .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
            discrete and continuous phase-type distributions of
            second order," in In Proc. of UK Performance 
            Evaluation Workshop, UKPEW, 2002"

    Examples
    --------
    For Matlab:
    
    >>> moms = [10,160,3500];
    >>> [a,A]=PH2From3Moments(moms);
    >>> a
       0.8702       0.1298
    >>> A
     -0.15576      0.15576
            0     -0.22659
    >>> MomentsFromPH(a,A,3)
          10         160        3500

