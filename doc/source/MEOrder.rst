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

    >>> a = [1./6, 1./6, 1./6, 1./6, 1./6, 1./6];
    >>> A = [-1., 0., 0., 0., 0., 0.; 0.5, -2., 1., 0., 0., 0.; 1., 0., -3., 1., 0., 0.; 1., 0., 1., -4., 1., 0.; 4., 0., 0., 0., -5., 0.; 5., 0., 0., 0., 0., -6.];
    >>> co = MEOrder(a,A,'cont');
    >>> disp(co);
         2
    >>> oo = MEOrder(a,A,'obs');
    >>> disp(oo);
         6
    >>> coo = MEOrder(a,A,'obscont');
    >>> disp(coo);
         2
    >>> mo = MEOrder(a,A,'moment');
    >>> disp(mo);
         2
    >>> a = [2./3, 1./3];
    >>> A = [-1., 1.; 0., -3.];
    >>> co = MEOrder(a,A,'cont');
    >>> disp(co);
         2
    >>> oo = MEOrder(a,A,'obs');
    >>> disp(oo);
         1
    >>> coo = MEOrder(a,A,'obscont');
    >>> disp(coo);
         1
    >>> mo = MEOrder(a,A,'moment');
    >>> disp(mo);
         1
    >>> b = [0.2, 0.3, 0.5];
    >>> B = [-1., 0., 0.; 0., -3., 1.; 0., -1., -3.];
    >>> [a,A] = MonocyclicPHFromME(b,B);
    >>> disp(a);
        0.0055089    0.0090301     0.016938     0.015216    0.0053543    0.0087356     0.052486      0.22657      0.66016
    >>> disp(A);
               -1            1            0            0            0            0            0            0            0
                0      -2.4226       2.4226            0            0            0            0            0            0
                0            0      -2.4226       2.4226            0            0            0            0            0
                0      0.26232            0      -2.4226       2.1603            0            0            0            0
                0            0            0            0      -4.2414       4.2414            0            0            0
                0            0            0            0            0      -4.2414       4.2414            0            0
                0            0            0            0            0            0      -4.2414       4.2414            0
                0            0            0            0            0            0            0      -4.2414       4.2414
                0            0            0            0            0            0            0            0      -4.2414
    >>> co = MEOrder(a,A,'cont');
    >>> disp(co);
         9
    >>> oo = MEOrder(a,A,'obs');
    >>> disp(oo);
         3
    >>> coo = MEOrder(a,A,'obscont');
    >>> disp(coo);
         3
    >>> mo = MEOrder(a,A,'moment');
    >>> disp(mo);
         3

    For Mathematica:

    >>> a = {1./6, 1./6, 1./6, 1./6, 1./6, 1./6};
    >>> A = {{-1., 0., 0., 0., 0., 0.},{0.5, -2., 1., 0., 0., 0.},{1., 0., -3., 1., 0., 0.},{1., 0., 1., -4., 1., 0.},{4., 0., 0., 0., -5., 0.},{5., 0., 0., 0., 0., -6.}};
    >>> co = MEOrder[a,A,"cont"];
    >>> Print[co];
    2
    >>> oo = MEOrder[a,A,"obs"];
    >>> Print[oo];
    6
    >>> coo = MEOrder[a,A,"obscont"];
    >>> Print[coo];
    2
    >>> mo = MEOrder[a,A,"moment"];
    >>> Print[mo];
    2
    >>> a = {2./3, 1./3};
    >>> A = {{-1., 1.},{0., -3.}};
    >>> co = MEOrder[a,A,"cont"];
    >>> Print[co];
    2
    >>> oo = MEOrder[a,A,"obs"];
    >>> Print[oo];
    1
    >>> coo = MEOrder[a,A,"obscont"];
    >>> Print[coo];
    1
    >>> mo = MEOrder[a,A,"moment"];
    >>> Print[mo];
    1
    >>> b = {0.2, 0.3, 0.5};
    >>> B = {{-1., 0., 0.},{0., -3., 1.},{0., -1., -3.}};
    >>> {a,A} = MonocyclicPHFromME[b,B];
    >>> Print[a];
    {0.00550893408977846, 0.00903007832853331, 0.016937512518639578, 0.015215980106503445, 0.005354337535618665, 0.008735592607040744, 0.05248568615571608, 0.22657249403204927, 0.6601593846261203}
    >>> Print[A];
    {{-1., 1., 0., 0., 0., 0., 0., 0., 0.},
     {0., -2.4226497308103743, 2.4226497308103743, 0., 0., 0., 0., 0., 0.},
     {0., 0., -2.4226497308103743, 2.4226497308103743, 0., 0., 0., 0., 0.},
     {0., 0.2623172489622428, 0., -2.4226497308103743, 2.1603324818481315, 0., 0., 0., 0.},
     {0., 0., 0., 0., -4.241399978863847, 4.241399978863847, 0., 0., 0.},
     {0., 0., 0., 0., 0., -4.241399978863847, 4.241399978863847, 0., 0.},
     {0., 0., 0., 0., 0., 0., -4.241399978863847, 4.241399978863847, 0.},
     {0., 0., 0., 0., 0., 0., 0., -4.241399978863847, 4.241399978863847},
     {0., 0., 0., 0., 0., 0., 0., 0., -4.241399978863847}}
    >>> co = MEOrder[a,A,"cont"];
    >>> Print[co];
    9
    >>> oo = MEOrder[a,A,"obs"];
    >>> Print[oo];
    3
    >>> coo = MEOrder[a,A,"obscont"];
    >>> Print[coo];
    3
    >>> mo = MEOrder[a,A,"moment"];
    >>> Print[mo];
    3

    For Python/Numpy:

    >>> a = ml.matrix([[1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])
    >>> A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])
    >>> co = MEOrder(a,A,"cont")
    >>> print(co)
    2
    >>> oo = MEOrder(a,A,"obs")
    >>> print(oo)
    6
    >>> coo = MEOrder(a,A,"obscont")
    >>> print(coo)
    2
    >>> mo = MEOrder(a,A,"moment")
    >>> print(mo)
    2
    >>> a = ml.matrix([[2./3, 1./3]])
    >>> A = ml.matrix([[-1., 1.],[0., -3.]])
    >>> co = MEOrder(a,A,"cont")
    >>> print(co)
    2
    >>> oo = MEOrder(a,A,"obs")
    >>> print(oo)
    1
    >>> coo = MEOrder(a,A,"obscont")
    >>> print(coo)
    1
    >>> mo = MEOrder(a,A,"moment")
    >>> print(mo)
    1
    >>> b = ml.matrix([[0.2, 0.3, 0.5]])
    >>> B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])
    >>> a,A = MonocyclicPHFromME(b,B)
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
    >>> co = MEOrder(a,A,"cont")
    >>> print(co)
    9
    >>> oo = MEOrder(a,A,"obs")
    >>> print(oo)
    3
    >>> coo = MEOrder(a,A,"obscont")
    >>> print(coo)
    3
    >>> mo = MEOrder(a,A,"moment")
    >>> print(mo)
    3

