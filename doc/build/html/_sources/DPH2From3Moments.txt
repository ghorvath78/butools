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
    ========
    For Matlab:

    >>> a = [0.9, 0.1];
    >>> A = [0.2, 0.61; 0.58, 0.41];
    >>> moms = MomentsFromDPH(a,A);
    >>> disp(moms);
           10.305       215.13       6764.2
    >>> [b,B] = DPH2From3Moments(moms);
    >>> disp(b);
          0.43249      0.56751
    >>> disp(B);
             0.61         0.39
          0.69692            0
    >>> phmoms = MomentsFromDPH(b,B,3);
    >>> disp(phmoms);
           10.305       215.13       6764.2

    For Mathematica:

    >>> a = {0.9, 0.1};
    >>> A = {{0.2, 0.61},{0.58, 0.41}};
    >>> moms = MomentsFromDPH[a,A];
    >>> Print[moms];
    {10.304568527918775, 215.1328300136563, 6764.166152521255}
    >>> {b,B} = DPH2From3Moments[moms];
    >>> Print[b];
    {0.43248730964467125, 0.5675126903553286}
    >>> Print[B];
    {{0.6100000000000014, 0.38999999999999857},
     {0.6969230769230763, 0}}
    >>> phmoms = MomentsFromDPH[b,B,3];
    >>> Print[phmoms];
    {10.304568527918788, 215.13283001365693, 6764.166152521286}

    For Python/Numpy:

    >>> a = ml.matrix([[0.9, 0.1]])
    >>> A = ml.matrix([[0.2, 0.61],[0.58, 0.41]])
    >>> moms = MomentsFromDPH(a,A)
    >>> print(moms)
    [10.304568527918775, 215.1328300136563, 6764.1661525212548]
    >>> b,B = DPH2From3Moments(moms)
    >>> print(b)
    [[ 0.43249  0.56751]]
    >>> print(B)
    [[ 0.61     0.39   ]
     [ 0.69692  0.     ]]
    >>> phmoms = MomentsFromDPH(b,B,3)
    >>> print(phmoms)
    [10.304568527918779, 215.1328300136565, 6764.1661525212639]

