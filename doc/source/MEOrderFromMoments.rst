butools.ph.MEOrderFromMoments
=============================

.. currentmodule:: butools.ph

.. np:function:: MEOrderFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`order = MEOrderFromMoments(moms, prec)`
        * - Mathematica:
          - :code:`order = MEOrderFromMoments[moms, prec]`
        * - Python/Numpy:
          - :code:`order = MEOrderFromMoments(moms, prec)`

    Returns the order of ME distribution that can realize
    the given moments.
    
    Parameters
    ----------
    moms : list of doubles
        The list of moments
    prec : double, optional
        Precision used to detect if the determinant of the
        Hankel matrix is zero. The default value is 1e-12.
    
    Returns
    -------
    order : int
        The order of ME distribution that can realize the 
        given moments

    References
    ----------
    .. [1]  L. Bodrog, A. Horvath, M. Telek, "Moment 
            characterization of matrix exponential and Markovian
            arrival processes," Annals of Operations Research, 
            vol. 160, pp. 51-68, 2008.

    Examples
    --------
    For Matlab:
    
    >>> a=[0.1 0.9 0];
    >>> A=[-6.2 2 0; 2 -9 1; 1 0 -3];
    >>> moms=MomentsFromME(a,A);
    >>> mo = MEOrderFromMoments(moms)
         3
     
    >>> b = [0.2, 0.3, 0.5];
    >>> B = [-1,0,0;0,-3,2;0,-2,-3];
    >>> [a,A] = MonocyclicPHFromME(b,B);
    >>> length(a)
        27
    >>> size(A)
        27    27
    >>> moms=MomentsFromME(a,A);
    >>> mo = MEOrderFromMoments(moms)
         3


