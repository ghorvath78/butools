butools.map.MarginalDistributionFromRAP
=======================================

.. currentmodule:: butools.map

.. np:function:: MarginalDistributionFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromRAP(H0, H1, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromRAP[H0, H1, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromRAP(H0, H1, precision)`

    Returns the phase type distributed marginal distribution
    of a rational arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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
    
    >>> H0=[-2 0 0; 0 -3 1; 0 -1 -2];
    >>> H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];
    >>> [a,A]=MarginalDistributionFromRAP(H0,H1);
    >>> a
          0.44444      0.44444      0.11111
    >>> A
        -2     0     0
         0    -3     1
         0    -1    -2
    >>> MomentsFromME(a,A)
          0.44444      0.38095      0.48299      0.82216       1.7944
    >>> x = (0:0.01:1);
    >>> plot(x,PdfFromME(a,A,x));

