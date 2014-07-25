butools.mam.FluidFundamentalMatrices
====================================

.. currentmodule:: butools.mam

.. np:function:: FluidFundamentalMatrices

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`M = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method)`
        * - Mathematica:
          - :code:`M = FluidFundamentalMatrices[xxx]`
        * - Python/Numpy:
          - :code:`M = FluidFundamentalMatrices(xxx)`

    Returns the fundamental matrices corresponding to the
    given canonical Markov fluid model. Matrices Psi, K and
    U are returned depending on the "matrices" parameter.
    
    Parameters
    ----------
    Fpp : matrix, shape (Np,Np)
        The matrix of transition rates between states 
        having positive fluid rates
    Fpm : matrix, shape (Np,Nm)
        The matrix of transition rates where the source
        state has a positive, the destination has a 
        negative fluid rate associated.
    Fpm : matrix, shape (Nm,Np)
        The matrix of transition rates where the source
        state has a negative, the destination has a 
        positive fluid rate associated.
    Fpp : matrix, shape (Nm,Nm)
        The matrix of transition rates between states 
        having negative fluid rates
    matrices : string
        Specifies which matrices are required. 'P' means 
        that only matrix Psi is needed. 'UK' means that
        matrices U and K are needed. Any combinations of
        'P', 'K' and 'U' are allowed, in any order.
    precision : double, optional
        The matrices are computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "ADDA", "SDA"}, optional
        The method used to solve the algebraic Riccati
        equation (CR: cyclic reduction, ADDA: alternating-
        directional doubling algorithm, SDA: structured
        doubling algorithm). The default is "CR".
    
    Returns
    -------
    M : list of matrices
        The list of calculated matrices in the order as
        requested in the 'matrices' parameter. If a single
        matrix is requested, the return value will be the
        matrix itself instead of a list.
    
    Examples
    --------
    For Matlab:
    
    >>> Fpp=[-5 1; 2 -3];
    >>> Fpm=[2 1 1; 1 0 0];
    >>> Fmm=[-8 4 1; 2 -12 3; 2 0 -2];
    >>> Fmp=[3 0; 2 5; 0 0];
    >>> M = FluidFundamentalMatrices(Fpp,Fpm,Fmp,Fmm,'PKU');
    >>> Psi=M{1}
          0.33722      0.16517      0.49761
           0.3318      0.12995      0.53825
    >>> K=M{2}
           -3.658       1.8258
           3.2553      -2.3502
    >>> U=M{3}
          -6.9883       4.4955       2.4928
           4.3334       -11.02       6.6865
                2            0           -2


