butools.dph.DPH2From3Moments
============================

.. currentmodule:: butools.dph

.. np:function:: DPH2From3Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = DPH2From3Moments(moms, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = DPH2From3Moments[moms, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = DPH2From3Moments(moms, prec)`

    Returns an order-2 discrete phase-type distribution 
    which has the same 3 moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,2)
      The initial probability vector of the DPH(2)
    A : matrix, shape (2,2)
      Transition probability matrix of the DPH(2)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a DPH(2).

    This procedure first calls 'MGFromMoments', then transforms
    it to DPH(2) by 'CanonicalFromDPH2'.

    Examples
    --------
    For Matlab:
    
    >>> moms = [10.305, 215.13, 6764.2];
    >>> [a,A]=DPH2From3Moments(moms)
    >>> a
      0.43249      0.56751
    >>> A
         0.61         0.39
      0.69692            0
    >>> MomentsFromDPH(a,A)
       10.305       215.13       6764.2

