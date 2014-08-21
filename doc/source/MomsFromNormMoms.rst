butools.moments.MomsFromNormMoms
================================

.. currentmodule:: butools.moments

.. np:function:: MomsFromNormMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`m = MomsFromNormMoms(nm)`
        * - Mathematica:
          - :code:`m = MomsFromNormMoms[nm]`
        * - Python/Numpy:
          - :code:`m = MomsFromNormMoms(nm)`
    
    Returns the raw moments given the normalized moments.
    
    The raw moments are: :math:`m_i=E(\mathcal{X}^i)`
    
    The normalized moments are: :math:`\displaystyle n_i=\frac{m_i}{m_{i-1} m_1}`
       
    Parameters
    ----------
    nm : vector of doubles
        The list of normalized moments (starting with the first
        moment)
        
    Returns
    -------
    m : vector of doubles
        The list of raw moments

    Examples
    --------
    For Matlab:
    
    >>> nm=NormMomsFromMoms([1.2 5 38 495 9215])
    [1.2 3.4722 6.3333 10.855 15.513]
    >>> m=MomsFromNormMoms(nm)
    [1.2 5 38 495 9215]
    
    For Mathematica:
    
    >>> nm=NormMomsFromMoms[{1.2, 5, 8}];
    >>> m=MomsFromNormMoms[nm]
    {1.2, 5, 8}
    
    For Python/Numpy:

    >>> nm=NormMomsFromMoms([1.2, 5, 8])
    >>> print(nm)
    [1.2, 3.4722222222222228, 6.333333333333333, 10.855263157894738, 15.513468013468014]
    >>> m=MomsFromNormMoms(nm)
    >>> print(m)
    [1.2, 5.000000000000001, 38.00000000000001, 495.00000000000017, 9215.000000000004]

