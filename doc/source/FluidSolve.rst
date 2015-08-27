butools.mam.FluidSolve
======================

.. currentmodule:: butools.mam

.. np:function:: FluidSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm, prec)`
        * - Mathematica:
          - :code:`{mass0, ini, K, clo} = FluidSolve [Fpp, Fpm, Fmp, Fmm, prec]`
        * - Python/Numpy:
          - :code:`mass0, ini, K, clo = FluidSolve (Fpp, Fpm, Fmp, Fmm, prec)`

    Returns the parameters of the matrix-exponentially 
    distributed stationary distribution of a canonical 
    Markovian fluid model.
    
    The canonical Markov fluid model is defined by the 
    matrix blocks of the generator of the background Markov
    chain partitioned according to the sign of the 
    associated fluid rates (i.e., there are "+" and "-" states).   

    Using the returned 4 parameters the stationary
    solution can be obtained as follows.
    
    The probability that the fluid level is zero while 
    being in different states of the background process
    is given by vector mass0.
    
    The density that the fluid level is x while being in
    different states of the background process is

    .. math::
        \pi(x)=ini\cdot e^{K x}\cdot clo.    
    
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
    precision : double, optional
        Numerical precision for computing the fundamental
        matrix. The default value is 1e-14
    
    Returns
    -------
    mass0 : matrix, shape (1,Np+Nm)
        The stationary probability vector of zero level
    ini : matrix, shape (1,Np)
        The initial vector of the stationary density
    K : matrix, shape (Np,Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np,Np+Nm)
        The closing matrix of the stationary density
    
    Examples
    ========
    For Matlab:

    >>> Fpp = [-5., 1.; 2., -3.];
    >>> Fpm = [2., 1., 1.; 1., 0., 0.];
    >>> Fmm = [-8., 4., 1.; 2., -12., 3.; 2., 0., -2.];
    >>> Fmp = [3., 0.; 2., 5.; 0., 0.];
    >>> x = 0.7;
    >>> [mass0, ini, K, clo] = FluidSolve(Fpp, Fpm, Fmp, Fmm);
    Final Residual Error for Psi:    1.1657e-15
    >>> disp(mass0);
         0.037514     0.015303     0.097921
    >>> disp(ini);
          0.14315     0.076517
    >>> disp(K);
           -3.658       1.8258
           3.2553      -2.3502
    >>> disp(clo);
                1            0      0.33722      0.16517      0.49761
                0            1       0.3318      0.12995      0.53825
    >>> pdfAtX = ini*expm(K*x)*clo;
    >>> disp(pdfAtX);
         0.074009     0.070933     0.048493     0.021442     0.075007

    For Mathematica:

    >>> Fpp = {{-5., 1.},{2., -3.}};
    >>> Fpm = {{2., 1., 1.},{1., 0., 0.}};
    >>> Fmm = {{-8., 4., 1.},{2., -12., 3.},{2., 0., -2.}};
    >>> Fmp = {{3., 0.},{2., 5.},{0., 0.}};
    >>> x = 0.7;
    >>> {mass0, ini, K, clo} = FluidSolve[Fpp, Fpm, Fmp, Fmm];
    "Final Residual Error for Psi: "1.1657341758564144*^-15
    >>> Print[mass0];
    {0.03751362697958451, 0.0153034356482914, 0.09792058971336806}
    >>> Print[ini];
    {0.14314775223533632, 0.076517178241457}
    >>> Print[K];
    {{-3.65799640319986, 1.8258294108775632},
     {3.255293764043593, -2.350237730252557}}
    >>> Print[clo];
    {{1, 0, 0.33722394414970486, 0.16516588217551262, 0.4976101736747833},
     {0, 1, 0.3317962853815385, 0.12995245394948857, 0.5382512606689742}}
    >>> pdfAtX = ini.MatrixExponential[K*x].clo;
    >>> Print[pdfAtX];
    {0.14314775223533632, 0.076517178241457} . MatrixExponential[{{-2.560597482239902, 1.2780805876142942}, {2.278705634830515, -1.6451664111767899}}] . {{1, 0, 0.33722394414970486, 0.16516588217551262, 0.4976101736747833}, {0, 1, 0.3317962853815385, 0.12995245394948857, 0.5382512606689742}}

    For Python/Numpy:

    >>> Fpp = ml.matrix([[-5., 1.],[2., -3.]])
    >>> Fpm = ml.matrix([[2., 1., 1.],[1., 0., 0.]])
    >>> Fmm = ml.matrix([[-8., 4., 1.],[2., -12., 3.],[2., 0., -2.]])
    >>> Fmp = ml.matrix([[3., 0.],[2., 5.],[0., 0.]])
    >>> x = 0.7
    >>> mass0, ini, K, clo = FluidSolve(Fpp, Fpm, Fmp, Fmm)
    Final Residual Error for G:  1.7208456881689926e-15
    >>> print(mass0)
    [[ 0.03751  0.0153   0.09792]]
    >>> print(ini)
    [[ 0.14315  0.07652]]
    >>> print(K)
    [[-3.658    1.82583]
     [ 3.25529 -2.35024]]
    >>> print(clo)
    [[ 1.       0.       0.33722  0.16517  0.49761]
     [ 0.       1.       0.3318   0.12995  0.53825]]
    >>> pdfAtX = ini*la.expm(K*x)*clo
    >>> print(pdfAtX)
    [[ 0.07401  0.07093  0.04849  0.02144  0.07501]]

