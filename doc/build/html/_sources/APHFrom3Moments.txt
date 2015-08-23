butools.ph.APHFrom3Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: APHFrom3Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = APHFrom3Moments(moms, maxSize)`
        * - Mathematica:
          - :code:`{alpha, A} = APHFrom3Moments[moms, maxSize]`
        * - Python/Numpy:
          - :code:`alpha, A = APHFrom3Moments(moms, maxSize)`

    Returns an acyclic PH which has the same 3 moments as
    given. If detects the order and the structure 
    automatically to match the given moments.

    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    maxSize : int, optional
      The maximal size of the resulting APH. The default value
      is 100.

    Returns
    -------
    alpha : vector, shape (1,M)
      The initial probability vector of the APH
    A : matrix, shape (M,M)
      Transient generator matrix of the APH
    
    Raises an error if the moments are not feasible with an
    APH of size "maxSize".
    
    References
    ----------
    .. [1] A. Bobbio, A. Horvath, M. Telek, "Matching three 
           moments with minimal acyclic phase type 
           distributions," Stochastic models, pp. 303-326, 2005.

    Examples
    ========
    For Matlab:

    >>> moms = [10.0, 125.0, 8400.0];
    >>> [a, A] = APHFrom3Moments(moms);
    >>> disp(a);
       1.3212e-05      0.99999            0            0            0            0
    >>> disp(A);
       -0.0022936    0.0022936            0            0            0            0
                0     -0.50029      0.50029            0            0            0
                0            0     -0.50029      0.50029            0            0
                0            0            0     -0.50029      0.50029            0
                0            0            0            0     -0.50029      0.50029
                0            0            0            0            0     -0.50029
    >>> phmoms = MomentsFromPH(a, A, 3);
    >>> disp(phmoms);
               10          125         8400
    >>> moms = [10.0, 525.0, 31400.0];
    >>> [a, A] = APHFrom3Moments(moms);
    >>> disp(a);
      Columns 1 through 6
          0.21179            0            0            0            0            0
      Columns 7 through 8
                0      0.78821
    >>> disp(A);
      Columns 1 through 6
         -0.15079      0.15079            0            0            0            0
                0     -0.15079      0.15079            0            0            0
                0            0     -0.15079      0.15079            0            0
                0            0            0     -0.15079      0.15079            0
                0            0            0            0     -0.15079      0.15079
                0            0            0            0            0     -0.15079
                0            0            0            0            0            0
                0            0            0            0            0            0
      Columns 7 through 8
                0            0
                0            0
                0            0
                0            0
                0            0
          0.15079            0
         -0.15079      0.15079
                0      -5.9502
    >>> phmoms = MomentsFromPH(a, A, 3);
    CheckMERepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!
    >>> disp(phmoms);
               10          525        31400

    For Mathematica:

    
    For Python/Numpy:

    >>> moms = [10.0, 125.0, 8400.0]
    >>> a, A = APHFrom3Moments(moms)
    >>> print(a)
    [[  1.32117e-05   9.99987e-01   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]]
    >>> print(A)
    [[-0.00229  0.00229  0.       0.       0.       0.     ]
     [ 0.      -0.50029  0.50029  0.       0.       0.     ]
     [ 0.       0.      -0.50029  0.50029  0.       0.     ]
     [ 0.       0.       0.      -0.50029  0.50029  0.     ]
     [ 0.       0.       0.       0.      -0.50029  0.50029]
     [ 0.       0.       0.       0.       0.      -0.50029]]
    >>> phmoms = MomentsFromPH(a, A, 3)
    >>> print(phmoms)
    [9.9999999999999964, 124.99999999999994, 8400.0000000001164]
    >>> moms = [10.0, 525.0, 31400.0]
    >>> a, A = APHFrom3Moments(moms)
    >>> print(a)
    [[ 0.21179  0.       0.       0.       0.       0.       0.       0.78821]]
    >>> print(A)
    [[-0.15079  0.15079  0.       0.       0.       0.       0.       0.     ]
     [ 0.      -0.15079  0.15079  0.       0.       0.       0.       0.     ]
     [ 0.       0.      -0.15079  0.15079  0.       0.       0.       0.     ]
     [ 0.       0.       0.      -0.15079  0.15079  0.       0.       0.     ]
     [ 0.       0.       0.       0.      -0.15079  0.15079  0.       0.     ]
     [ 0.       0.       0.       0.       0.      -0.15079  0.15079  0.     ]
     [ 0.       0.       0.       0.       0.       0.      -0.15079  0.15079]
     [ 0.       0.       0.       0.       0.       0.       0.      -5.9502 ]]
    >>> phmoms = MomentsFromPH(a, A, 3)
    CheckMERepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!
    >>> print(phmoms)
    [10.0, 525.0, 31399.999999999993]

