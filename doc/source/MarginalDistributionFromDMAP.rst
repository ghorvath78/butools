butools.dmap.MarginalDistributionFromDMAP
=========================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalDistributionFromDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromDMAP(D0, D1, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromDMAP[D0, D1, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromDMAP(D0, D1, precision)`

    Returns the discrete phase type distributed marginal 
    distribution of a discrete Markovian arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase 
        type distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the discrete phase type 
        distributed marginal distribution    

    Examples
    --------
    For Matlab:
    
    >>> D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12];
    >>> D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27];
    >>> [a,A]=MarginalDistributionFromDMAP(D0,D1);
    >>> a
          0.24388      0.40412       0.1941       0.1579
    >>> A
                0         0.02            0            0
                0         0.17          0.2         0.14
             0.16         0.17         0.02         0.18
                0            0            0         0.12
    >>> MomentsFromDPH(a,A)
           1.4955       2.9542       7.8852       27.282       116.17       587.04         3437
    >>> x = (1:10);
    >>> bar(x,PmfFromDPH(a,A,x));

