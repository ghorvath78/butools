butools.dmap.MarginalDistributionFromDMRAP
==========================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalDistributionFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromDMRAP(H, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromDMRAP[H, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromDMRAP(H, precision)`

    Returns the matrix geometrically distributed marginal 
    distribution of a discrete marked rational arrival 
    process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix geometrically
        distributed marginal distribution
    A : matrix, shape (M,M)
        The matrix parameter of the matrix geometrically
        distributed marginal distribution    

    Examples
    --------
    For Matlab:
    
    >>> H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16];
    >>> H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17];
    >>> H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0];
    >>> [a,A]=MarginalDistributionFromDMRAP({H0,H1,H2});
    >>> a
          0.22615      0.35424      0.41962
    >>> A
             0.15          0.2         0.18
            -0.23         0.17         0.22
             0.19         0.15         0.16
    >>> MomentsFromMG(a,A)
           1.5948       3.4185       9.9595       37.742       177.13
    >>> x = (1:10);
    >>> bar(x,PmfFromMG(a,A,x));

    For Python/Numpy:

    >>> H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    >>> [a,A]=MarginalDistributionFromDMRAP((H0,H1,H2))
    >>> print(a)
    [[ 0.22614581  0.35423782  0.41961638]]
    >>> print(A)
    [[ 0.15  0.2   0.18]
     [-0.23  0.17  0.22]
     [ 0.19  0.15  0.16]]
    >>> print(MomentsFromMG(a,A))
    [1.5947682612697014, 3.4185320639932879, 9.9594855925243557, 37.742369501645818, 177.12501064870213]
    >>> x = np.linspace(1,10,10)
    >>> plt.bar(x,PmfFromMG(a,A,x))

