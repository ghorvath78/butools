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
    ========
    For Matlab:

    >>> M = [1.2, 5., 38., 495., 9215.];
    >>> nmoms = NormMomsFromMoms(M);
    >>> disp(nmoms);
              1.2       3.4722       6.3333       10.855       15.513
    >>> moms = MomsFromNormMoms(nmoms);
    >>> disp(moms);
              1.2            5           38          495         9215
    >>> err = norm(moms-M);
    >>> disp(err);
         0

    For Mathematica:

    >>> M = {1.2, 5., 38., 495., 9215.};
    >>> nmoms = NormMomsFromMoms[M];
    >>> Print[nmoms];
    {1.2, 3.4722222222222223, 6.333333333333333, 10.855263157894736, 15.513468013468012}
    >>> moms = MomsFromNormMoms[nmoms];
    >>> Print[moms];
    {1.2, 5., 37.99999999999999, 494.99999999999983, 9214.999999999996}
    >>> err = Norm[moms-M];
    >>> Print[err];
    3.641980348156267*^-12

    For Python/Numpy:

    >>> M = [1.2, 5., 38., 495., 9215.]
    >>> nmoms = NormMomsFromMoms(M)
    >>> print(nmoms)
    [1.2, 3.4722222222222228, 6.333333333333333, 10.855263157894738, 15.513468013468014]
    >>> moms = MomsFromNormMoms(nmoms)
    >>> print(moms)
    [1.2, 5.000000000000001, 38.00000000000001, 495.00000000000017, 9215.000000000004]
    >>> err = la.norm(np.array(moms)-np.array(M))
    >>> print(err)
    3.641980456457359e-12

