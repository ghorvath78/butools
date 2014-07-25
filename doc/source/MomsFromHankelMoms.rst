butools.moments.MomsFromHankelMoms
==================================

.. currentmodule:: butools.moments

.. np:function:: MomsFromHankelMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`hm = MomsFromHankelMoms(m)`
        * - Mathematica:
          - :code:`hm = MomsFromHankelMoms[m]`
        * - Python/Numpy:
          - :code:`hm = MomsFromHankelMoms(m)`
    
    Returns the raw moments given the Hankel moments.
    
    The raw moments are: :math:`m_i=E(\mathcal{X}^i)`
    
    The ith Hankel moment is the determinant of matrix 
    :math:`\Delta_{i/2}`, if i is even, 
    and it is the determinant of :math:`\Delta^{(1)}_{(i+1)/2}`, 
    if i is odd. For the definition of matrices :math:`\Delta`
    and :math:`\Delta^{(1)}` see [1]_.
       
    Parameters
    ----------
    hm : vector of doubles
        The list of Hankel moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem

    Examples
    --------
    For Matlab:
    
    >>> hm=HankelMomsFromMoms([1.3 2.4 6.03 20.5 89.5 474.9])
    [1.3 0.71 2.079 1.9973 13.841 44.916]
    >>> m=MomsFromHankelMoms(hm)
    [1.3 2.4 6.03 20.5 89.5 474.9]
    
    For Mathematica:
    
    >>> hm=HankelMomsFromMoms[{1.2, 5, 38}];
    >>> m=MomsFromHankelMoms[hm]
    {1.2, 5, 38}
    
    For Python/Numpy:

    >>> hm=HankelMomsFromMoms([1.2, 5, 38]);
    >>> m=MomsFromHankelMoms(hm)
    [1.2, 5, 38]

