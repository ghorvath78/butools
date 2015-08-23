butools.ph.MEOrder
==================

.. currentmodule:: butools.ph

.. np:function:: MEOrder

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`order = MEOrder(alpha, A, kind, prec)`
        * - Mathematica:
          - :code:`order = MEOrder[alpha, A, kind, prec]`
        * - Python/Numpy:
          - :code:`order = MEOrder(alpha, A, kind, prec)`

    Returns the order of the ME distribution (which is not 
    necessarily equal to the size of the representation).

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    kind : {'obs', 'cont', 'obscont', 'moment'}, optional
        Determines which order is computed. Possibilities: 
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
        controllability order,
        'moment': moment order (which is the default).
    prec : double, optional
        Precision used to detect if the determinant of the 
        Hankel matrix is zero (in case of kind="moment" only),
        or the tolerance for the rank calculation. The
        default value is 1e-10.

    Returns
    -------
    order : int
        The order of ME distribution

    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation
            of rational arrival processes." Madrid Conference on
            Qeueuing theory (MCQT), June 2010.

    Examples
    ========
    For Matlab:

    >>> a = [1.0/6,1.0/6,1.0/6,1.0/6,1.0/6,1.0/6];
    >>> A = [-1., 0., 0., 0., 0., 0.; 0.5, -2., 1., 0., 0., 0.; 1., 0., -3., 1., 0., 0.; 1., 0., 1., -4., 1., 0.; 4., 0., 0., 0., -5., 0.; 5., 0., 0., 0., 0., -6.];
    >>> co = MEOrder(a, A, 'cont');
    >>> disp(co);
         2
    >>> oo = MEOrder(a, A, 'obs');
    >>> disp(oo);
         6
    >>> coo = MEOrder(a, A, 'obscont');
    >>> disp(coo);
         2
    >>> mo = MEOrder(a, A, 'moment');
    >>> disp(mo);
         2
    >>> a = [2.0/3,1.0/3];
    >>> A = [-1., 1.; 0., -3.];
    >>> co = MEOrder(a, A, 'cont');
    >>> disp(co);
         2
    >>> oo = MEOrder(a, A, 'obs');
    >>> disp(oo);
         1
    >>> coo = MEOrder(a, A, 'obscont');
    >>> disp(coo);
         1
    >>> mo = MEOrder(a, A, 'moment');
    >>> disp(mo);
         1
    >>> b = [0.2,0.3,0.5];
    >>> B = [-1., 0., 0.; 0., -3., 1.; 0., -1., -3.];
    >>> [a, A] = MonocyclicPHFromME(b, B);
    >>> disp(a);
      Columns 1 through 6
        0.0055089    0.0090301     0.016938     0.015216    0.0053543    0.0087356
      Columns 7 through 9
         0.052486      0.22657      0.66016
    >>> disp(A);
      Columns 1 through 6
               -1            1            0            0            0            0
                0      -2.4226       2.4226            0            0            0
                0            0      -2.4226       2.4226            0            0
                0      0.26232            0      -2.4226       2.1603            0
                0            0            0            0      -4.2414       4.2414
                0            0            0            0            0      -4.2414
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
      Columns 7 through 9
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
           4.2414            0            0
          -4.2414       4.2414            0
                0      -4.2414       4.2414
                0            0      -4.2414
    >>> co = MEOrder(a, A, 'cont');
    >>> disp(co);
         9
    >>> oo = MEOrder(a, A, 'obs');
    >>> disp(oo);
         3
    >>> coo = MEOrder(a, A, 'obscont');
    >>> disp(coo);
         3
    >>> mo = MEOrder(a, A, 'moment');
    >>> disp(mo);
         3

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[1.0/6,1.0/6,1.0/6,1.0/6,1.0/6,1.0/6]])
    >>> A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])
    >>> co = MEOrder(a, A, "cont")
    >>> print(co)
    2
    >>> oo = MEOrder(a, A, "obs")
    >>> print(oo)
    6
    >>> coo = MEOrder(a, A, "obscont")
    >>> print(coo)
    2
    >>> mo = MEOrder(a, A, "moment")
    >>> print(mo)
    2
    >>> a = ml.matrix([[2.0/3,1.0/3]])
    >>> A = ml.matrix([[-1., 1.],[0., -3.]])
    >>> co = MEOrder(a, A, "cont")
    >>> print(co)
    2
    >>> oo = MEOrder(a, A, "obs")
    >>> print(oo)
    1
    >>> coo = MEOrder(a, A, "obscont")
    >>> print(coo)
    1
    >>> mo = MEOrder(a, A, "moment")
    >>> print(mo)
    1
    >>> b = ml.matrix([[0.2,0.3,0.5]])
    >>> B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])
    >>> a, A = MonocyclicPHFromME(b, B)
    >>> print(a)
    [[ 0.00551  0.00903  0.01694  0.01522  0.00535  0.00874  0.05249  0.22657  0.66016]]
    >>> print(A)
    [[-1.       1.       0.       0.       0.       0.       0.       0.       0.     ]
     [ 0.      -2.42265  2.42265  0.       0.       0.       0.       0.       0.     ]
     [ 0.       0.      -2.42265  2.42265  0.       0.       0.       0.       0.     ]
     [ 0.       0.26232  0.      -2.42265  2.16033  0.       0.       0.       0.     ]
     [ 0.       0.       0.       0.      -4.2414   4.2414   0.       0.       0.     ]
     [ 0.       0.       0.       0.       0.      -4.2414   4.2414   0.       0.     ]
     [ 0.       0.       0.       0.       0.       0.      -4.2414   4.2414   0.     ]
     [ 0.       0.       0.       0.       0.       0.       0.      -4.2414   4.2414 ]
     [ 0.       0.       0.       0.       0.       0.       0.       0.      -4.2414 ]]
    >>> co = MEOrder(a, A, "cont")
    >>> print(co)
    9
    >>> oo = MEOrder(a, A, "obs")
    >>> print(oo)
    3
    >>> coo = MEOrder(a, A, "obscont")
    >>> print(coo)
    3
    >>> mo = MEOrder(a, A, "moment")
    >>> print(mo)
    3

