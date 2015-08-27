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
    ========
    For Matlab:

    >>> Fpp = [-5., 1.; 2., -3.];
    >>> Fpm = [2., 1., 1.; 1., 0., 0.];
    >>> Fmm = [-8., 4., 1.; 2., -12., 3.; 2., 0., -2.];
    >>> Fmp = [3., 0.; 2., 5.; 0., 0.];
    >>> [Psi, K, U] = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, 'PKU');
    Final Residual Error for Psi:    1.1657e-15
    >>> disp(Psi);
          0.33722      0.16517      0.49761
           0.3318      0.12995      0.53825
    >>> disp(K);
           -3.658       1.8258
           3.2553      -2.3502
    >>> disp(U);
          -6.9883       4.4955       2.4928
           4.3334       -11.02       6.6865
                2            0           -2

    For Mathematica:

    >>> Fpp = {{-5., 1.},{2., -3.}};
    >>> Fpm = {{2., 1., 1.},{1., 0., 0.}};
    >>> Fmm = {{-8., 4., 1.},{2., -12., 3.},{2., 0., -2.}};
    >>> Fmp = {{3., 0.},{2., 5.},{0., 0.}};
    >>> {Psi, K, U} = FluidFundamentalMatrices[Fpp, Fpm, Fmp, Fmm, "PKU"];
    "Final Residual Error for Psi: "1.1657341758564144*^-15
    >>> Print[Psi];
    {{0.33722394414970486, 0.16516588217551262, 0.4976101736747833},
     {0.3317962853815385, 0.12995245394948857, 0.5382512606689742}}
    >>> Print[K];
    {{-3.65799640319986, 1.8258294108775632},
     {3.255293764043593, -2.350237730252557}}
    >>> Print[U];
    {{-6.988328167550885, 4.4954976465265375, 2.4928305210243495},
     {4.333429315207102, -11.019905965901533, 6.686476650694438},
     {2., 0., -2.}}

    For Python/Numpy:

    >>> Fpp = ml.matrix([[-5., 1.],[2., -3.]])
    >>> Fpm = ml.matrix([[2., 1., 1.],[1., 0., 0.]])
    >>> Fmm = ml.matrix([[-8., 4., 1.],[2., -12., 3.],[2., 0., -2.]])
    >>> Fmp = ml.matrix([[3., 0.],[2., 5.],[0., 0.]])
    >>> Psi, K, U = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, "PKU")
    Final Residual Error for G:  1.7208456881689926e-15
    >>> print(Psi)
    [[ 0.33722  0.16517  0.49761]
     [ 0.3318   0.12995  0.53825]]
    >>> print(K)
    [[-3.658    1.82583]
     [ 3.25529 -2.35024]]
    >>> print(U)
    [[ -6.98833   4.4955    2.49283]
     [  4.33343 -11.01991   6.68648]
     [  2.        0.       -2.     ]]

