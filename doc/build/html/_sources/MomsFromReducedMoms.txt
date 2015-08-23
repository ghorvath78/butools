butools.moments.MomsFromReducedMoms
===================================

.. currentmodule:: butools.moments

.. np:function:: MomsFromReducedMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m = MomsFromReducedMoms(rm)`
        * - Mathematica:
          - :code:`m = MomsFromReducedMoms[rm]`
        * - Python/Numpy:
          - :code:`m = MomsFromReducedMoms(rm)`
    
    Returns the raw moments given the reduced moments.
    
    The raw moments are: :math:`m_i=E(\mathcal{X}^i)`
    
    The reduced moments are: :math:`\displaystyle r_i=\frac{m_i}{i!}`
       
    Parameters
    ----------
    rm : vector of doubles
        The list of reduced moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments

    Examples
    ========
    For Matlab:

    >>> M = [1.2, 5., 38., 495., 9215.];
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

    
    For Python/Numpy:

    >>> M = [1.2, 5., 38., 495., 9215.]
    >>> rmoms = ReducedMomsFromMoms(M)
    >>> print(rmoms)
    [1.2, 2.5, 6.333333333333333, 20.625, 76.79166666666667]
    >>> moms = MomsFromReducedMoms(rmoms)
    >>> print(moms)
    [1.2, 5.0, 38.0, 495.0, 9215.0]
    >>> err = la.norm(np.array(moms)-np.array(M))
    >>> print(err)
    0.0

