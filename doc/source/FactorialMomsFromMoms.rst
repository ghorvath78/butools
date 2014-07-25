butools.moments.FactorialMomsFromMoms
=====================================

.. currentmodule:: butools.moments

.. np:function:: FactorialMomsFromMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`fm = FactorialMomsFromMoms(m)`
        * - Mathematica:
          - :code:`fm = FactorialMomsFromMoms[m]`
        * - Python/Numpy:
          - :code:`fm = FactorialMomsFromMoms(m)`
    
    Returns the factorial moments given the raw moments.
    
    The raw moments are: :math:`m_i=E(\mathcal{X}^i)`
    
    The factorial moments are: :math:`f_i=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1))`
       
    Parameters
    ----------
    m : vector of doubles
        The list of raw moments (starting with the first
        moment)
        
    Returns
    -------
    fm : vector of doubles
        The list of factorial moments

    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    

    Examples
    --------
    For Matlab:
    
    >>> fm=FactorialMomsFromMoms([1.3 2.4 6.03 20.5 89.5 474.9])
    [1.3 1.1 1.43 2.92 6.75 19.75]
    >>> m=MomsFromFactorialMoms(fm)
    [1.3 2.4 6.03 20.5 89.5 474.9]

    For Mathematica:
    
    >>> fm=FactorialMomsFromMoms[{1.2, 5, 38}]
    
    For Python/Numpy:

    >>> fm=FactorialMomsFromMoms([1.2, 5, 38])

