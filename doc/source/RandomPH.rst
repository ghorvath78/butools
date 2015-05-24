butools.ph.RandomPH
===================

.. currentmodule:: butools.ph

.. np:function:: RandomPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = RandomPH(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = RandomPH[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = RandomPH(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random phase-type distribution with a given 
    order.

    Parameters
    ----------
    order : int
        The size of the phase-type distribution
    mean : double, optional
        The mean of the phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper PH 
        (that has an irreducible phase process and none of 
        its parameters is all-zero). The default value is 
        1000.
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type 
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type 
        distribution.

    Notes
    -----
    If the procedure fails, try to increase the 'maxTrials'
    parameter.   

    Examples
    ========
    For Matlab:

    >>> [a,A] = RandomPH(3,8,4);
    >>> disp(a);
         0     1     0
    >>> disp(A);
         -0.97262      0.56772      0.14713
          0.23072     -0.69732      0.46659
          0.54311            0     -0.71127

    For Mathematica:

    >>> {a,A} = RandomPH[3,8,4];
    >>> Print[a];
    {0.10568235872400128, 0.2894271458934268, 0.604890495382572}
    >>> Print[A];
    {{-0.18520978401771204, 0., 0.},
     {0., -0.1156138846071312, 0.},
     {0.2582356194770228, 0.08094833396062599, -0.37998447481873276}}

    For Python/Numpy:

    >>> a,A = RandomPH(3,8,4)
    >>> print(a)
    [[ 0.51378  0.48622  0.     ]]
    >>> print(A)
    [[-0.68768  0.14149  0.29254]
     [ 0.06547 -0.21547  0.     ]
     [ 0.       0.16149 -0.16149]]

