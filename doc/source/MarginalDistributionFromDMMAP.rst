butools.dmap.MarginalDistributionFromDMMAP
==========================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalDistributionFromDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromDMMAP(D, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromDMMAP[D, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromDMMAP(D, precision)`

    Returns the discrete phase type distributed marginal 
    distribution of a discrete marked Markovian arrival 
    process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
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

    >>> D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0];
    >>> D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09];
    >>> D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13];
    >>> D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01];
    >>> [a,A]=MarginalDistributionFromDMMAP({D0,D1,D2,D3});
    >>> a
          0.96166     0.030652    0.0076858
    >>> A
             0.34            0            0
             0.06         0.05         0.03
             0.11         0.13            0
    >>> MomentsFromDPH(a,A)
           1.5037       3.0278       8.4243       31.097       143.88
    >>> x = (1:10);
    >>> bar(x,PmfFromDPH(a,A,x));

    For Python/Numpy:

    >>> D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    >>> D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    >>> D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    >>> D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    >>> [a,A]=MarginalDistributionFromDMMAP((D0,D1,D2,D3))
    >>> print(a)
    [[ 0.96166245  0.03065179  0.00768576]]
    >>> print(A)
    [[ 0.34  0.    0.  ]
     [ 0.06  0.05  0.03]
     [ 0.11  0.13  0.  ]]
    >>> print(MomentsFromDPH(a,A))
    [1.503697331491185, 3.0278418573508938, 8.424305390832199, 31.097386717744087, 143.88041840101141]
    >>> x = np.linspace(1,10,10)
    >>> plt.bar(x,PmfFromDPH(a,A,x))

