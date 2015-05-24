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
    --------
    For Matlab:

    >>> [alpha,A]=RandomPH(3,8,4);
    >>> alpha
      0.79483     0.072975      0.13219
    >>> A
      -3.7518       3.7518            0
       19.001      -36.518       15.859
       6.2631            0      -6.2631
    >>> sum(A,2)
            0
      -1.6586
            0
    >>> mean = MomentsFromPH (alpha, A, 1)
            8

    For Python/Numpy:

    >>> alpha,A = RandomPH(3,8,4)
    >>> print(alpha)
    [[ 0.53928457  0.46071543  0.        ]]
    >>> print(A)
    [[-1.14323874  0.1175927   0.46955765]
     [ 0.59310796 -1.23831621  0.64520825]
     [ 0.          0.6843269  -0.6843269 ]]
    >>> print(np.sum(A,1))
    [[-0.55608839]
     [ 0.        ]
     [ 0.        ]]
    >>> print(MomentsFromPH (alpha, A, 1))
    [7.9999999999999911]

