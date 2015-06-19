butools.map.RandomMAP
=====================

.. currentmodule:: butools.map

.. np:function:: RandomMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = RandomMAP(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{D0, D1} = RandomMAP[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D0, D1 = RandomMAP(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random Markovian arrival process with given mean 
    value.

    Parameters
    ----------
    order : int
        The size of the MAP
    mean : double, optional
        The mean inter-arrival times of the MAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper MAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D0 : vector, shape (1,M)
        The D0 matrix of the MAP
    D1 : matrix, shape (M,M)
        The D1 matrix of the MAP

    Examples
    ========
    For Matlab:

    >>> [D0,D1] = RandomMAP(4,1.62,10);
    >>> disp(D0);
          -2.9064            0      0.42583      0.51808
                0      -2.6084            0      0.43424
          0.10604      0.58736      -1.0914            0
          0.28667            0      0.32664     -0.61331
    >>> disp(D1);
          0.64541       0.0731      0.51853      0.72543
          0.63998      0.55402      0.27188      0.70824
          0.10328      0.28404     0.010665            0
                0            0            0            0
    >>> m = MarginalMomentsFromMAP(D0,D1,1);
    >>> disp(m);
             1.62

    For Mathematica:

    >>> {D0,D1} = RandomMAP[4,1.62,10];
    >>> Print[D0];
    {{-2.0923812019623003, 0.13304510567104486, 0.31829701350927053, 0.1521668074411143},
     {0., -0.5741730780189052, 0., 0.},
     {0., 0.34151979439016145, -1.1427165394526269, 0.13142372139224592},
     {0.40438386153011857, 0.07922513291605328, 0., -1.6331535398111858}}
    >>> Print[D1];
    {{0.2271800914455872, 0.29319862534449953, 0.4661257840587526, 0.5023677744920312},
     {0., 0.5373250633033104, 0., 0.03684801471559468},
     {0., 0., 0.20387931620380573, 0.46589370746641395},
     {0.3566073647001836, 0.40462364617714686, 0., 0.38831353448768324}}
    >>> m = MarginalMomentsFromMAP[D0,D1,1][[1]];
    >>> Print[m];
    1.6199999999999999

    For Python/Numpy:

    >>> D0,D1 = RandomMAP(4,1.62,10)
    >>> print(D0)
    [[-2.43095  0.19245  0.48386  0.2269 ]
     [ 0.18365 -2.37753  0.52155  0.41858]
     [ 0.       0.54859 -1.52778  0.58001]
     [ 0.       0.19143  0.      -0.69681]]
    >>> print(D1)
    [[ 0.47978  0.63877  0.22953  0.17967]
     [ 0.44676  0.40001  0.       0.40698]
     [ 0.       0.       0.39918  0.     ]
     [ 0.       0.       0.       0.50538]]
    >>> m = MarginalMomentsFromMAP(D0,D1,1)[0]
    >>> print(m)
    1.62

