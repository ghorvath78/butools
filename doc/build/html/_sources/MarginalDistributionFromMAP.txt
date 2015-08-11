butools.map.MarginalDistributionFromMAP
=======================================

.. currentmodule:: butools.map

.. np:function:: MarginalDistributionFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromMAP(D0, D1, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromMAP[D0, D1, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromMAP(D0, D1, precision)`

    Returns the phase type distributed marginal distribution
    of a Markovian arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase type 
        distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the phase type distributed
        marginal distribution    

    Examples
    ========
    For Matlab:

    >>> D0 = [-0.17, 0, 0, 0.07; 0.01, -0.78, 0.03, 0.08; 0.22, 0.17, -1.1, 0.02; 0.04, 0.12, 0, -0.42];
    >>> D1 = [0, 0.06, 0, 0.04; 0.04, 0.19, 0.21, 0.22; 0.22, 0.13, 0.15, 0.19; 0.05, 0, 0.17, 0.04];
    >>> [a,A] = MarginalDistributionFromMAP(D0,D1);
    >>> disp(a);
          0.14438      0.23571      0.33794      0.28196
    >>> disp(A);
            -0.17            0            0         0.07
             0.01        -0.78         0.03         0.08
             0.22         0.17         -1.1         0.02
             0.04         0.12            0        -0.42

    For Mathematica:

    >>> D0 = {{-0.17, 0, 0, 0.07},{0.01, -0.78, 0.03, 0.08},{0.22, 0.17, -1.1, 0.02},{0.04, 0.12, 0, -0.42}};
    >>> D1 = {{0, 0.06, 0, 0.04},{0.04, 0.19, 0.21, 0.22},{0.22, 0.13, 0.15, 0.19},{0.05, 0, 0.17, 0.04}};
    >>> {a,A} = MarginalDistributionFromMAP[D0,D1];
    >>> Print[a];
    {0.1443829740435819, 0.23570990103732542, 0.33794226877724004, 0.2819648561418526}
    >>> Print[A];
    {{-0.17, 0, 0, 0.07},
     {0.01, -0.78, 0.03, 0.08},
     {0.22, 0.17, -1.1, 0.02},
     {0.04, 0.12, 0, -0.42}}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    >>> D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    >>> a,A = MarginalDistributionFromMAP(D0,D1)
    >>> print(a)
    [[ 0.14438  0.23571  0.33794  0.28196]]
    >>> print(A)
    [[-0.17  0.    0.    0.07]
     [ 0.01 -0.78  0.03  0.08]
     [ 0.22  0.17 -1.1   0.02]
     [ 0.04  0.12  0.   -0.42]]

