butools.mam.FluidFundamentalMatrices
====================================

.. currentmodule:: butools.mam

.. np:function:: FluidFundamentalMatrices

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`M = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method)`
        * - Mathematica:
          - :code:`M = FluidFundamentalMatrices[Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method]`
        * - Python/Numpy:
          - :code:`M = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method)`

    Returns the fundamental matrices corresponding to the
    given canonical Markov fluid model. Matrices Psi, K and
    U are returned depending on the "matrices" parameter.
    The canonical Markov fluid model is defined by the 
    matrix blocks of the generator of the background Markov
    chain partitioned according to the sign of the 
    associated fluid rates (i.e., there are "+" and "-" states).
    
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
        requested in the 'matrices' parameter.
    
    Notes
    -----
    Thanks to Benny Van Houdt for the implementation of the
    Riccati solvers.
    
    Examples
    --------
    For Matlab:
    
    >>> Fpp=[-5 1; 2 -3];
    >>> Fpm=[2 1 1; 1 0 0];
    >>> Fmm=[-8 4 1; 2 -12 3; 2 0 -2];
    >>> Fmp=[3 0; 2 5; 0 0];
    >>> [Psi, K, U] = FluidFundamentalMatrices(Fpp,Fpm,Fmp,Fmm,'PKU');
    >>> Psi
          0.33722      0.16517      0.49761
           0.3318      0.12995      0.53825
    >>> K
           -3.658       1.8258
           3.2553      -2.3502
    >>> U
          -6.9883       4.4955       2.4928
           4.3334       -11.02       6.6865
                2            0           -2

    For Python/Numpy:

    >>> Fpp=ml.matrix([[-5, 1],[2, -3]])
    >>> Fpm=ml.matrix([[2, 1, 1],[1, 0, 0]])
    >>> Fmm=ml.matrix([[-8, 4, 1],[2, -12, 3],[2, 0, -2]])
    >>> Fmp=ml.matrix([[3, 0],[2, 5],[0, 0]])
    >>> Psi, K, U = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, "PKU")
    >>> print(Psi)
    [[ 0.33722394  0.16516588  0.49761017]
     [ 0.33179629  0.12995245  0.53825126]]
    >>> print(K)
    [[-3.6579964   1.82582941]
     [ 3.25529376 -2.35023773]]    
    >>> print(U)
    [[ -6.98832817   4.49549765   2.49283052]
     [  4.33342932 -11.01990597   6.68647665]
     [  2.           0.          -2.        ]]

