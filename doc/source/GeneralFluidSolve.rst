butools.mam.GeneralFluidSolve
=============================

.. currentmodule:: butools.mam

.. np:function:: GeneralFluidSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[mass0, ini, K, clo] = GeneralFluidSolve (Q, R, Q0)`
        * - Mathematica:
          - :code:`{mass0, ini, K, clo} = GeneralFluidSolve [Q, R, Q0]`
        * - Python/Numpy:
          - :code:`mass0, ini, K, clo = GeneralFluidSolve (Q, R, Q0)`

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
    
    >>> Q=[-6 1 3 2 0 0; 6 -10 2 0 2 0; 3 7 -12 0 0 2; 5 0 0 -9 1 3; 0 5 0 6 -13 2; 0 0 5 3 7 -15];
    >>> R=[2 0 0 0 0 0; 0 -4 0 0 0 0; 0 0 -12 0 0 0; 0 0 0 6 0 0; 0 0 0 0 0 0; 0 0 0 0 0 -8];
    >>> [mass0,ini,K,clo]=GeneralFluidSolve(Q,R);
    >>> mass0
                0     0.082246     0.069492            0     0.023812     0.020724
    >>> ini
          0.70195      0.20505
    >>> K
          -2.4698       1.1349
            1.295      -1.1686
    >>> clo
              0.5     0.061087     0.054574            0      0.01618     0.012595
                0     0.055389     0.043116      0.16667     0.038913     0.032631
    >>> emptyprob = sum(mass0)
          0.19628
    >>> densityat1 = ini*exp(K*1)*clo
      0.40402      0.17384      0.14099      0.37455      0.10052      0.08351


