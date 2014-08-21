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
    --------
    For Matlab:
    
    >>> rm=ReducedMomsFromMoms([1.2 5 38 495 9215])
    [1.2 2.5 6.3333 20.625 76.792]
    >>> m=MomsFromReducedMoms(rm)
    [1.2 5 38 495 9215]
    
    For Mathematica:
    
    >>> rm=ReducedMomsFromMoms[{1.2, 5, 8}];
    >>> m=MomsFromReducedMoms[rm]
    {1.2, 5, 8}
    
    For Python/Numpy:

    >>> rm=ReducedMomsFromMoms([1.2, 5, 38, 495, 9215])
    >>> print(rm)
    [1.2, 2.5, 6.333333333333333, 20.625, 76.79166666666667]
    >>> m=MomsFromReducedMoms(rm)
    >>> print(m)
    [1.2, 5.0, 38.0, 495.0, 9215.0]

