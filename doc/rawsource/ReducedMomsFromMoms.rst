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
    --------
    For Matlab:
    
    >>> rm=ReducedMomsFromMoms([1.2 5 38 495 9215])
    [1.2 2.5 6.3333 20.625 76.792]
    >>> m=MomsFromReducedMoms(rm)
    [1.2 5 38 495 9215]
    
    For Mathematica:
    
    >>> rm=ReducedMomsFromMoms[{1.2, 5, 38, 495, 9215}]
    {1.2, 5/2, 19/3, 165/8, 1843/24}
    >>> m=MomsFromReducedMoms[rm]
    {1.2, 5, 38, 495, 9215}
    
    For Python/Numpy:

    >>> rm=ReducedMomsFromMoms([1.2, 5, 38, 495, 9215])
    >>> print(rm)
    [1.2, 2.5, 6.333333333333333, 20.625, 76.79166666666667]
    >>> m=MomsFromReducedMoms(rm)
    >>> print(m)
    [1.2, 5.0, 38.0, 495.0, 9215.0]

