butools.ph.AcyclicPHFromME
==========================

.. currentmodule:: butools.ph

.. np:function:: AcyclicPHFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = AcyclicPHFromME(alpha, A, maxSize, precision)`
        * - Mathematica:
          - :code:`{beta, B} = AcyclicPHFromME[alpha, A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`beta, B = AcyclicPHFromME(alpha, A, maxSize, precision)`

    Transforms an arbitrary matrix-exponential representation
    to an acyclic phase-type representation. (see [1]_).
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    maxSize : int, optional
        The maximum number of phases for the result.
        The default value is 100.
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the Markovian 
        acyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        acyclic representation
    
    Notes
    -----
    Raises an error if no Markovian acyclic representation
    has been found.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
            phase-type distributions," Stoch. Models 15, 759-778 
            (1999)

    Examples
    ========
    For Matlab:

    >>> a = [-0.4, 1.4, 0.];
    >>> A = [-4., 1., 1.; 0., -2., 1.; 1., 0., -8.];
    >>> [b,B] = AcyclicPHFromME(a,A);
    >>> disp(b);
          0.55273       0.3741     0.073173
    >>> disp(B);
          -1.9145       1.9145            0
                0      -3.8858       3.8858
                0            0      -8.1997
    >>> ma = MomentsFromME(a,A,5);
    >>> disp(ma);
          0.64918      0.73131       1.1825       2.5062       6.5898
    >>> mb = MomentsFromME(b,B,5);
    >>> disp(mb);
          0.64918      0.73131       1.1825       2.5062       6.5898

    For Mathematica:

    >>> a = {-0.4, 1.4, 0.};
    >>> A = {{-4., 1., 1.},{0., -2., 1.},{1., 0., -8.}};
    >>> {b,B} = AcyclicPHFromME[a,A];
    >>> Print[b];
    {0.5527262576934738, 0.3741003779708655, 0.07317336433566028}
    >>> Print[B];
    {{-1.914468283493477, 1.914468283493477, 0},
     {0, -3.8858267357749496, 3.8858267357749496},
     {0, 0, -8.199704980731571}}
    >>> ma = MomentsFromME[a,A,5];
    >>> Print[ma];
    {0.6491803278688524, 0.7313087879602256, 1.1825377454500594, 2.50620910640242, 6.589761685446926}
    >>> mb = MomentsFromME[b,B,5];
    >>> Print[mb];
    {0.6491803278688524, 0.7313087879602258, 1.18253774545006, 2.506209106402422, 6.5897616854469305}

    For Python/Numpy:

    >>> a = ml.matrix([[-0.4, 1.4, 0.]])
    >>> A = ml.matrix([[-4., 1., 1.],[0., -2., 1.],[1., 0., -8.]])
    >>> b,B = AcyclicPHFromME(a,A)
    >>> print(b)
    [[ 0.55273  0.3741   0.07317]]
    >>> print(B)
    [[-1.91447  1.91447  0.     ]
     [ 0.      -3.88583  3.88583]
     [ 0.       0.      -8.1997 ]]
    >>> ma = MomentsFromME(a,A,5)
    >>> print(ma)
    [0.64918032786885238, 0.73130878796022558, 1.1825377454500594, 2.5062091064024199, 6.5897616854469261]
    >>> mb = MomentsFromME(b,B,5)
    >>> print(mb)
    [0.64918032649380331, 0.73130878479925954, 1.182537737569604, 2.5062090835980859, 6.5897616090589066]

