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
         0.098961     0.027218     0.022577     0.045517      0.01383     0.011405

    For Python/Numpy:

    >>> Q=ml.matrix([[-6, 1, 3, 2, 0, 0],[6, -10, 2, 0, 2, 0],[3, 7, -12, 0, 0, 2],[5, 0, 0, -9, 1, 3],[0, 5, 0, 6, -13, 2],[0, 0, 5, 3, 7, -15]])
    >>> R=ml.matrix([[2, 0, 0, 0, 0, 0],[0, -4, 0, 0, 0, 0],[0, 0, -12, 0, 0, 0],[0, 0, 0, 6, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, -8]])
    >>> mass0,ini,K,clo=GeneralFluidSolve(Q,R)
    >>> print(mass0)
    [[ 0.          0.08224613  0.0694924   0.          0.02381248  0.02072428]]
    >>> print(ini)
    [[ 0.70195397  0.20504772]]
    >>> print(K)
    [[-2.46975215  1.1348626 ]
     [ 1.29501774 -1.16863056]]
    >>> print(clo)
    [[ 0.5         0.06108742  0.05457444  0.          0.01617979  0.01259463]
     [ 0.          0.05538938  0.04311605  0.16666667  0.03891262  0.03263124]]
    >>> emptyprob = np.sum(mass0)
    >>> print(emptyprob)
    0.196275290213
    >>> densityat1 = ini*la.expm(K*1)*clo
    >>> print(densityat1)
    [[ 0.09896148  0.0272177   0.02257673  0.04551746  0.01382957  0.01140451]]

