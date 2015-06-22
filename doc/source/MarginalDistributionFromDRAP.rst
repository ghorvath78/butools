butools.dmap.MarginalDistributionFromDRAP
=========================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalDistributionFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromDRAP(H0, H1, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromDRAP[H0, H1, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromDRAP(H0, H1, precision)`

    Returns the matrix geometrically distributed marginal 
    distribution of a discrete rational arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
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
    ========
    For Matlab:

    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1, -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> [a,A] = MarginalDistributionFromDRAP(H0,H1);
    >>> disp(a);
         0.021493      0.71253      0.26598
    >>> disp(A);
                0            0         0.13
                0          0.6         0.18
             0.31         0.26         0.02

    For Mathematica:

    >>> H0 = {{0, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1, -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> {a,A} = MarginalDistributionFromDRAP[H0,H1];
    >>> Print[a];
    {0.021493050580312315, 0.7125271919655068, 0.26597975745418073}
    >>> Print[A];
    {{0, 0, 0.13},
     {0, 0.6, 0.18},
     {0.31, 0.26, 0.02}}

    For Python/Numpy:

    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> a,A = MarginalDistributionFromDRAP(H0,H1)
    >>> print(a)
    [[ 0.02149  0.71253  0.26598]]
    >>> print(A)
    [[ 0.    0.    0.13]
     [ 0.    0.6   0.18]
     [ 0.31  0.26  0.02]]

