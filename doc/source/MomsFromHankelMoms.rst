butools.moments.MomsFromHankelMoms
==================================

.. currentmodule:: butools.moments

.. np:function:: MomsFromHankelMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m = MomsFromHankelMoms(hm)`
        * - Mathematica:
          - :code:`m = MomsFromHankelMoms[hm]`
        * - Python/Numpy:
          - :code:`m = MomsFromHankelMoms(hm)`
    
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
    
    >>> hm=HankelMomsFromMoms[{1.3, 2.4, 6.03, 20.5, 89.5, 474.9}];
    {1.3, 0.71, 2.079, 1.9973, 13.8413, 44.9157}
    >>> m=MomsFromHankelMoms[hm]
    {1.3, 2.4, 6.03, 20.5, 89.5, 474.9}
    
    For Python/Numpy:

    >>> hm=HankelMomsFromMoms([1.3, 2.4, 6.03, 20.5, 89.5, 474.9])
    >>> print(hm)
    [1.3, 0.7099999999999995, 2.079000000000002, 1.997299999999989, 13.841272999999994, 44.91574881000027]
    >>> m=MomsFromHankelMoms(hm)
    >>> print(m)
    [1.3, 2.3999999999999995, 6.0299999999999985, 20.499999999999996, 89.500000000000057, 474.90000000000032]

