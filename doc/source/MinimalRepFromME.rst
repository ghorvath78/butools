butools.ph.MinimalRepFromME
===========================

.. currentmodule:: butools.ph

.. np:function:: MinimalRepFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = MinimalRepFromME(alpha, A, how, precision)`
        * - Mathematica:
          - :code:`{beta, B} = MinimalRepFromME[alpha, A, how, precision]`
        * - Python/Numpy:
          - :code:`beta, B = MinimalRepFromME(alpha, A, how, precision)`

    Returns the minimal representation of the given ME 
    distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    how : {"obs", "cont", "obscont", "moment"}, optional        
        Determines how the representation is minimized. 
        Possibilities:
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
        controllability order,
        'moment': moment order (which is the default).
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.

    Returns
    -------
    beta : vector, shape (1,N)
        The initial vector of the minimal representation
    B : matrix, shape (N,N)
        The matrix parameter of the minimal representation

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
    >>> [b, B] = MinimalRepFromME(a, A, 'cont');
    >>> disp(b);
                1   1.3878e-16
    >>> disp(B);
          -1.4011      0.48448
          0.49585      -1.5989
    >>> [b, B] = MinimalRepFromME(a, A, 'obs');
    >>> disp(b);
          0.16667      0.16667      0.16667      0.16667      0.16667      0.16667
    >>> disp(B);
               -1            0            0            0            0            0
              0.5           -2            1            0            0            0
                1            0           -3            1            0            0
                1            0            1           -4            1            0
                4            0            0            0           -5            0
                5            0            0            0            0           -6
    >>> [b, B] = MinimalRepFromME(a, A, 'obscont');
    >>> disp(b);
                1   1.3878e-16
    >>> disp(B);
          -1.4011      0.48448
          0.49585      -1.5989
    >>> [b, B] = MinimalRepFromME(a, A, 'moment');
    >>> disp(b);
              0.5          0.5
    >>> disp(B);
            -2.52       1.6467
            -0.48        -0.48
    >>> a = [2.0/3,1.0/3];
    >>> A = [-1., 1.; 0., -3.];
    >>> [b, B] = MinimalRepFromME(a, A, 'cont');
    >>> disp(b);
          0.66667      0.33333
    >>> disp(B);
        -1     1
         0    -3
    >>> [b, B] = MinimalRepFromME(a, A, 'obs');
    >>> disp(b);
                1
    >>> disp(B);
        -1
    >>> [b, B] = MinimalRepFromME(a, A, 'obscont');
    >>> disp(b);
                1
    >>> disp(B);
        -1
    >>> [b, B] = MinimalRepFromME(a, A, 'moment');
    >>> disp(b);
         1
    >>> disp(B);
        -1
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
    >>> [b, B] = MinimalRepFromME(a, A, 'cont');
    >>> disp(b);
      Columns 1 through 6
        0.0055089    0.0090301     0.016938     0.015216    0.0053543    0.0087356
      Columns 7 through 9
         0.052486      0.22657      0.66016
    >>> disp(B);
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
    >>> [b, B] = MinimalRepFromME(a, A, 'obs');
    >>> disp(b);
                1   2.0817e-17  -5.5511e-17
    >>> disp(B);
          -2.8362     0.036222  -4.4409e-16
           -16.61      -3.3369       16.042
           1.1643    -0.051724     -0.82688
    >>> Cm = SimilarityMatrix(B, A);
    >>> err1 = norm(B*Cm-Cm*A);
    >>> err2 = norm(b*Cm-a);
    >>> disp(max(err1, err2));
        9.334e-15
    >>> [b, B] = MinimalRepFromME(a, A, 'obscont');
    >>> disp(b);
                1   2.0817e-17  -5.5511e-17
    >>> disp(B);
          -2.8362     0.036222  -4.4409e-16
           -16.61      -3.3369       16.042
           1.1643    -0.051724     -0.82688
    >>> Cm = SimilarityMatrix(B, A);
    >>> err1 = norm(B*Cm-Cm*A);
    >>> err2 = norm(b*Cm-a);
    >>> disp(max(err1, err2));
        9.334e-15
    >>> [b, B] = MinimalRepFromME(a, A, 'moment');
    >>> disp(b);
          0.33333      0.33333      0.33333
    >>> disp(B);
          -2.1905       1.9222      -3.3698
          -1.0769      -2.3906      0.83162
         -0.51037       0.8033      -2.4189
    >>> Cm = SimilarityMatrix(B, A);
    >>> err1 = norm(B*Cm-Cm*A);
    >>> err2 = norm(b*Cm-a);
    >>> disp(max(err1, err2));
       5.5343e-15

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[1.0/6,1.0/6,1.0/6,1.0/6,1.0/6,1.0/6]])
    >>> A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])
    >>> b, B = MinimalRepFromME(a, A, "cont")
    >>> print(b)
    [[  1.00000e+00   2.08167e-16]]
    >>> print(B)
    [[-1.40115  0.48448]
     [ 0.49585 -1.59885]]
    >>> b, B = MinimalRepFromME(a, A, "obs")
    >>> print(b)
    [[ 0.16667  0.16667  0.16667  0.16667  0.16667  0.16667]]
    >>> print(B)
    [[-1.   0.   0.   0.   0.   0. ]
     [ 0.5 -2.   1.   0.   0.   0. ]
     [ 1.   0.  -3.   1.   0.   0. ]
     [ 1.   0.   1.  -4.   1.   0. ]
     [ 4.   0.   0.   0.  -5.   0. ]
     [ 5.   0.   0.   0.   0.  -6. ]]
    >>> b, B = MinimalRepFromME(a, A, "obscont")
    >>> print(b)
    [[  1.00000e+00   2.08167e-16]]
    >>> print(B)
    [[-1.40115  0.48448]
     [ 0.49585 -1.59885]]
    >>> b, B = MinimalRepFromME(a, A, "moment")
    >>> print(b)
    [[ 0.5  0.5]]
    >>> print(B)
    [[-2.52     1.64667]
     [-0.48    -0.48   ]]
    >>> a = ml.matrix([[2.0/3,1.0/3]])
    >>> A = ml.matrix([[-1., 1.],[0., -3.]])
    >>> b, B = MinimalRepFromME(a, A, "cont")
    >>> print(b)
    [[ 0.66667  0.33333]]
    >>> print(B)
    [[-1.  1.]
     [ 0. -3.]]
    >>> b, B = MinimalRepFromME(a, A, "obs")
    >>> print(b)
    [[ 1.]]
    >>> print(B)
    [[-1.]]
    >>> b, B = MinimalRepFromME(a, A, "obscont")
    >>> print(b)
    [[ 1.]]
    >>> print(B)
    [[-1.]]
    >>> b, B = MinimalRepFromME(a, A, "moment")
    >>> print(b)
    [[ 1.]]
    >>> print(B)
    [[-1.]]
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
    >>> b, B = MinimalRepFromME(a, A, "cont")
    >>> print(b)
    [[ 0.00551  0.00903  0.01694  0.01522  0.00535  0.00874  0.05249  0.22657  0.66016]]
    >>> print(B)
    [[-1.       1.       0.       0.       0.       0.       0.       0.       0.     ]
     [ 0.      -2.42265  2.42265  0.       0.       0.       0.       0.       0.     ]
     [ 0.       0.      -2.42265  2.42265  0.       0.       0.       0.       0.     ]
     [ 0.       0.26232  0.      -2.42265  2.16033  0.       0.       0.       0.     ]
     [ 0.       0.       0.       0.      -4.2414   4.2414   0.       0.       0.     ]
     [ 0.       0.       0.       0.       0.      -4.2414   4.2414   0.       0.     ]
     [ 0.       0.       0.       0.       0.       0.      -4.2414   4.2414   0.     ]
     [ 0.       0.       0.       0.       0.       0.       0.      -4.2414   4.2414 ]
     [ 0.       0.       0.       0.       0.       0.       0.       0.      -4.2414 ]]
    >>> b, B = MinimalRepFromME(a, A, "obs")
    >>> print(b)
    [[  1.00000e+00   6.93889e-18   5.55112e-17]]
    >>> print(B)
    [[ -2.83622e+00   3.62221e-02  -2.22045e-16]
     [ -1.66096e+01  -3.33690e+00   1.60422e+01]
     [  1.16435e+00  -5.17236e-02  -8.26876e-01]]
    >>> Cm = SimilarityMatrix(B, A)
    >>> err1 = la.norm(B*Cm-Cm*A)
    >>> err2 = la.norm(b*Cm-a)
    >>> print(np.max(err1, err2))
    7.41724165831e-15
    >>> b, B = MinimalRepFromME(a, A, "obscont")
    >>> print(b)
    [[  1.00000e+00   6.93889e-18   5.55112e-17]]
    >>> print(B)
    [[ -2.83622e+00   3.62221e-02  -2.22045e-16]
     [ -1.66096e+01  -3.33690e+00   1.60422e+01]
     [  1.16435e+00  -5.17236e-02  -8.26876e-01]]
    >>> Cm = SimilarityMatrix(B, A)
    >>> err1 = la.norm(B*Cm-Cm*A)
    >>> err2 = la.norm(b*Cm-a)
    >>> print(np.max(err1, err2))
    7.41724165831e-15
    >>> b, B = MinimalRepFromME(a, A, "moment")
    >>> print(b)
    [[ 0.33333  0.33333  0.33333]]
    >>> print(B)
    [[-2.19048  1.92219 -3.36981]
     [-1.07693 -2.3906   0.83162]
     [-0.51037  0.8033  -2.41893]]
    >>> Cm = SimilarityMatrix(B, A)
    >>> err1 = la.norm(B*Cm-Cm*A)
    >>> err2 = la.norm(b*Cm-a)
    >>> print(np.max(err1, err2))
    2.3964518895e-15

