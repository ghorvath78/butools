butools.mam.FluidStationaryDistr
================================

.. currentmodule:: butools.mam

.. np:function:: FluidStationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = FluidStationaryDistr (mass0, ini, K, clo, x)`
        * - Mathematica:
          - :code:`pi = FluidStationaryDistr [mass0, ini, K, clo, x]`
        * - Python/Numpy:
          - :code:`pi = FluidStationaryDistr (mass0, ini, K, clo, x)`

    Returns the stationary distribution of a Markovian 
    fluid model at the given points.
    
    Parameters
    ----------
    mass0 : matrix, shape (1,Np+Nm)
        The stationary probability vector of zero level
    ini : matrix, shape (1,Np)
        The initial vector of the stationary density
    K : matrix, shape (Np,Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np,Np+Nm)
        The closing matrix of the stationary density
    x : vector, length (K)
        The distribution function is computed at these 
        points.
    
    Returns
    -------
    pi : matrix, shape (K,Nm+Np)
        The ith row of pi is the probability that the fluid
        level is less than or equal to x(i), while being in
        different states of the background process.
    
    Examples
    --------
    For Matlab:
    
    >>> Q=[-6 1 3 2 0 0; 6 -10 2 0 2 0; 3 7 -12 0 0 2; 5 0 0 -9 1 3; 0 5 0 6 -13 2; 0 0 5 3 7 -15];
    >>> R=[2 0 0 0 0 0; 0 -4 0 0 0 0; 0 0 -12 0 0 0; 0 0 0 6 0 0; 0 0 0 0 0 0; 0 0 0 0 0 -8];
    >>> [mass0,ini,K,clo]=GeneralFluidSolve(Q,R);
    >>> pi = FluidStationaryDistr (mass0,ini,K,clo, [0, 10, 100]);
    >>> pi(1,:)
                0     0.082246     0.069492            0     0.023812     0.020724
    >>> pi(2,:)
          0.37951      0.17891      0.15007      0.15134     0.071428     0.059915
    >>> pi(3,:)
          0.38328      0.18002      0.15099      0.15331     0.072009     0.060395
    >>> CTMCSolve(Q)
          0.38328      0.18002      0.15099      0.15331     0.072009     0.060395

    For Python/Numpy:

    >>> Q=ml.matrix([[-6, 1, 3, 2, 0, 0],[6, -10, 2, 0, 2, 0],[3, 7, -12, 0, 0, 2],[5, 0, 0, -9, 1, 3],[0, 5, 0, 6, -13, 2],[0, 0, 5, 3, 7, -15]])
    >>> R=ml.matrix([[2, 0, 0, 0, 0, 0],[0, -4, 0, 0, 0, 0],[0, 0, -12, 0, 0, 0],[0, 0, 0, 6, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, -8]])
    >>> mass0,ini,K,clo=GeneralFluidSolve(Q,R)
    >>> pi = FluidStationaryDistr (mass0,ini,K,clo, [0, 10, 100])
    >>> print(pi[0,:])
    [[ 0.          0.08224613  0.0694924   0.          0.02381248  0.02072428]]
    >>> print(pi[1,:])
    [[ 0.37950645  0.17890948  0.15006733  0.15134435  0.07142838  0.05991509]]
    >>> print(pi[2,:])
    [[ 0.38327526  0.18002323  0.15098722  0.1533101   0.07200929  0.06039489]]
    >>> print(CTMCSolve(Q))
    [[ 0.38327526  0.18002323  0.15098722  0.1533101   0.07200929  0.06039489]]

