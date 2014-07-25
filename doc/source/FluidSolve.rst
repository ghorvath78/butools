butools.mam.FluidSolve
======================

.. currentmodule:: butools.mam

.. np:function:: FluidSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm)`
        * - Mathematica:
          - :code:`{mass0, ini, K, clo} = FluidSolve [Fpp, Fpm, Fmp, Fmm]`
        * - Python/Numpy:
          - :code:`mass0, ini, K, clo = FluidSolve (Fpp, Fpm, Fmp, Fmm)`

    Returns the parameters of the matrix-exponentially 
    distributed stationary distribution of a canonical 
    Markovian fluid model

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
    --------
    For Matlab:
    
    >>> Fpp=[-5 1; 2 -3];
    >>> Fpm=[2 1 1; 1 0 0];
    >>> Fmm=[-8 4 1; 2 -12 3; 2 0 -2];
    >>> Fmp=[3 0; 2 5; 0 0];
    >>> [mass0,ini,K,clo]=FluidSolve(Fpp,Fpm,Fmp,Fmm);
    >>> mass0
         0.037514     0.015303     0.097921
    >>> ini
          0.14315     0.076517
    >>> K
           -3.658       1.8258
           3.2553      -2.3502
    >>> clo
                1            0      0.33722      0.16517      0.49761
                0            1       0.3318      0.12995      0.53825
    >>> emptyprob = sum(mass0)
          0.15074    
    >>> densityat1 = ini*exp(K*1)*clo
           1.9876      0.89595      0.96753      0.44471       1.4713
  

