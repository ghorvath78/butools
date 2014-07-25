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
    --------
    For Matlab:
    
    >>> H0=[-5 0.28 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
    >>> H1=[-0.08 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
    >>> H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]
    >>> [a,A]=MarginalDistributionFromMRAP({H0,H1,H2});
    >>> a
          0.17159      0.21695      0.27936       0.3321
    >>> A
               -5         0.28          0.9            1
                1           -8          0.9          0.1
              0.9          0.1           -4            1
                1            2            3           -9
    >>> MomentsFromME(a,A)
          0.33951      0.24583      0.27424      0.41206      0.77677       1.7594       4.6515
    >>> x = (0:0.01:1);
    >>> plot(x,PdfFromME(a,A,x));

