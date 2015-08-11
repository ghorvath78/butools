butools.map.MarginalDistributionFromMMAP
========================================

.. currentmodule:: butools.map

.. np:function:: MarginalDistributionFromMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromMMAP(D, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromMMAP[D, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromMMAP(D, precision)`

    Returns the phase type distributed marginal distribution
    of a marked Markovian arrival process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
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

    >>> D0 = [-1.78, 0.29; 0.07, -0.92];
    >>> D1 = [0.15, 0.49; 0.23, 0.36];
    >>> D2 = [0.11, 0.2; 0.01, 0];
    >>> D3 = [0.14, 0.4; 0.11, 0.14];
    >>> [a,A] = MarginalDistributionFromMMAP({D0,D1,D2,D3});
    >>> disp(a);
          0.36191      0.63809
    >>> disp(A);
            -1.78         0.29
             0.07        -0.92

    For Mathematica:

    >>> D0 = {{-1.78, 0.29},{0.07, -0.92}};
    >>> D1 = {{0.15, 0.49},{0.23, 0.36}};
    >>> D2 = {{0.11, 0.2},{0.01, 0}};
    >>> D3 = {{0.14, 0.4},{0.11, 0.14}};
    >>> {a,A} = MarginalDistributionFromMMAP[{D0,D1,D2,D3}];
    >>> Print[a];
    {0.36190793862575055, 0.6380920613742495}
    >>> Print[A];
    {{-1.78, 0.29},
     {0.07, -0.92}}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    >>> D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    >>> D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
    >>> D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    >>> a,A = MarginalDistributionFromMMAP([D0,D1,D2,D3])
    >>> print(a)
    [[ 0.36191  0.63809]]
    >>> print(A)
    [[-1.78  0.29]
     [ 0.07 -0.92]]

