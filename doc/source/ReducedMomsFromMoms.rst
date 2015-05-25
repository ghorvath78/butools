butools.moments.ReducedMomsFromMoms
===================================

.. currentmodule:: butools.moments

.. np:function:: ReducedMomsFromMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`rm = ReducedMomsFromMoms(m)`
        * - Mathematica:
          - :code:`rm = ReducedMomsFromMoms[m]`
        * - Python/Numpy:
          - :code:`rm = ReducedMomsFromMoms(m)`
    
    Returns the reduced moments given the raw moments.
    
    The raw moments are: :math:`m_i=E(\mathcal{X}^i)`
    
    The reduced moments are: :math:`\displaystyle r_i=\frac{m_i}{i!}`
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    rm : vector of doubles
        The list of reduced moments

    Examples
    ========
    For Matlab:

    >>> M = [1.2, 5, 38, 495, 9215];
    >>> rmoms = ReducedMomsFromMoms(M);
    >>> disp(rmoms);
              1.2          2.5       6.3333       20.625       76.792
    >>> moms = MomsFromReducedMoms(rmoms);
    >>> disp(moms);
              1.2            5           38          495         9215
    >>> err = norm(moms-M);
    >>> disp(err);
         0

    For Mathematica:

    >>> M = {1.2, 5, 38, 495, 9215};
    >>> rmoms = ReducedMomsFromMoms[M];
    >>> Print[rmoms];
    {1.2, 5/2, 19/3, 165/8, 1843/24}
    >>> moms = MomsFromReducedMoms[rmoms];
    >>> Print[moms];
    {1.2, 5, 38, 495, 9215}
    >>> err = Norm[moms-M];
    >>> Print[err];
    0.

    For Python/Numpy:

    >>> M = [1.2, 5, 38, 495, 9215]
    >>> rmoms = ReducedMomsFromMoms(M)
    >>> print(rmoms)
    [1.2, 2.5, 6.333333333333333, 20.625, 76.79166666666667]
    >>> moms = MomsFromReducedMoms(rmoms)
    >>> print(moms)
    [1.2, 5.0, 38.0, 495.0, 9215.0]
    >>> err = la.norm(np.array(moms)-np.array(M))
    >>> print(err)
    0.0

