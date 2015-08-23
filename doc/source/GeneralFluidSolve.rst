butools.mam.GeneralFluidSolve
=============================

.. currentmodule:: butools.mam

.. np:function:: GeneralFluidSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[mass0, ini, K, clo] = GeneralFluidSolve (Q, R, Q0, prec)`
        * - Mathematica:
          - :code:`{mass0, ini, K, clo} = GeneralFluidSolve [Q, R, Q0, prec]`
        * - Python/Numpy:
          - :code:`mass0, ini, K, clo = GeneralFluidSolve (Q, R, Q0, prec)`

    Returns the parameters of the matrix-exponentially 
    distributed stationary distribution of a general 
    Markovian fluid model, where the fluid rates associated
    with the states of the background process can be
    arbitrary (zero is allowed as well).

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
    Q : matrix, shape (N,N)
        The generator of the background Markov chain
    R : diagonal matrix, shape (N,N)
        The diagonal matrix of the fluid rates associated
        with the different states of the background process
    Q0 : matrix, shape (N,N), optional
        The generator of the background Markov chain at 
        level 0. If not provided, or empty, then Q0=Q is 
        assumed. The default value is empty.
    precision : double, optional
        Numerical precision for computing the fundamental
        matrix and for checking. The default value is 1e-14
    
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

    >>> Q = [-6., 1., 3., 2., 0., 0.; 6., -10., 2., 0., 2., 0.; 3., 7., -12., 0., 0., 2.; 5., 0., 0., -9., 1., 3.; 0., 5., 0., 6., -13., 2.; 0., 0., 5., 3., 7., -15.];
    >>> R = [2., 0., 0., 0., 0., 0.; 0., -4., 0., 0., 0., 0.; 0., 0., -12., 0., 0., 0.; 0., 0., 0., 6., 0., 0.; 0., 0., 0., 0., 0., 0.; 0., 0., 0., 0., 0., -8.];
    >>> x = 0.7;
    >>> [mass0, ini, K, clo] = GeneralFluidSolve(Q, R);
    Final Residual Error for Psi:    7.6328e-16
    >>> disp(mass0);
                0     0.082246     0.069492            0     0.023812     0.020724
    >>> disp(ini);
          0.70195      0.20505
    >>> disp(K);
          -2.4698       1.1349
            1.295      -1.1686
    >>> disp(clo);
              0.5     0.061087     0.054574            0      0.01618     0.012595
                0     0.055389     0.043116      0.16667     0.038913     0.032631
    >>> pdfAtX = ini*expm(K*x)*clo;
    >>> disp(pdfAtX);
          0.12566     0.031849     0.026557     0.049637     0.015655     0.012884

    For Mathematica:

    
    For Python/Numpy:

    >>> Q = ml.matrix([[-6., 1., 3., 2., 0., 0.],[6., -10., 2., 0., 2., 0.],[3., 7., -12., 0., 0., 2.],[5., 0., 0., -9., 1., 3.],[0., 5., 0., 6., -13., 2.],[0., 0., 5., 3., 7., -15.]])
    >>> R = ml.matrix([[2., 0., 0., 0., 0., 0.],[0., -4., 0., 0., 0., 0.],[0., 0., -12., 0., 0., 0.],[0., 0., 0., 6., 0., 0.],[0., 0., 0., 0., 0., 0.],[0., 0., 0., 0., 0., -8.]])
    >>> x = 0.7
    >>> mass0, ini, K, clo = GeneralFluidSolve(Q, R)
    Final Residual Error for G:  6.661338147750939e-16
    >>> print(mass0)
    [[ 0.       0.08225  0.06949  0.       0.02381  0.02072]]
    >>> print(ini)
    [[ 0.70195  0.20505]]
    >>> print(K)
    [[-2.46975  1.13486]
     [ 1.29502 -1.16863]]
    >>> print(clo)
    [[ 0.5      0.06109  0.05457  0.       0.01618  0.01259]
     [ 0.       0.05539  0.04312  0.16667  0.03891  0.03263]]
    >>> pdfAtX = ini*la.expm(K*x)*clo
    >>> print(pdfAtX)
    [[ 0.12566  0.03185  0.02656  0.04964  0.01566  0.01288]]

