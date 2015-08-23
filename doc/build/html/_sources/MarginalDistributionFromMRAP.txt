butools.map.MarginalDistributionFromMRAP
========================================

.. currentmodule:: butools.map

.. np:function:: MarginalDistributionFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromMRAP(H, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromMRAP[H, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromMRAP(H, precision)`

    Returns the phase type distributed marginal distribution
    of a marked rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix exponentially
        distributed marginal
    A : matrix, shape (M,M)
        The matrix parameter of the matrix exponentially
        distributed marginal    

    Examples
    ========
    For Matlab:

    >>> x = 0.18;
    >>> H0 = [-5., 0.1+x, 0.9, 1.; 1., -8., 0.9, 0.1; 0.9, 0.1, -4., 1.; 1., 2., 3., -9.];
    >>> H1 = [0.1-x, 0.7, 0.1, 0.1; 0.1, 1., 1.8, 0.1; 0.1, 0.1, 0.1, 0.7; 0.7, 0.1, 0.1, 0.1];
    >>> H2 = [0.1, 0.1, 0.1, 1.7; 1.8, 0.1, 1., 0.1; 0.1, 0.1, 0.7, 0.1; 0.1, 1., 0.1, 0.8];
    >>> [a, A] = MarginalDistributionFromMRAP({H0, H1, H2});
    >>> disp(a);
          0.17159      0.21695      0.27936       0.3321
    >>> disp(A);
               -5         0.28          0.9            1
                1           -8          0.9          0.1
              0.9          0.1           -4            1
                1            2            3           -9

    For Mathematica:

    
    For Python/Numpy:

    >>> x = 0.18
    >>> H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
    >>> H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
    >>> a, A = MarginalDistributionFromMRAP([H0, H1, H2])
    >>> print(a)
    [[ 0.17159  0.21695  0.27936  0.3321 ]]
    >>> print(A)
    [[-5.    0.28  0.9   1.  ]
     [ 1.   -8.    0.9   0.1 ]
     [ 0.9   0.1  -4.    1.  ]
     [ 1.    2.    3.   -9.  ]]

