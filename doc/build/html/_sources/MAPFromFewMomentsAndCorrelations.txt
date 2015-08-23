butools.map.MAPFromFewMomentsAndCorrelations
============================================

.. currentmodule:: butools.map

.. np:function:: MAPFromFewMomentsAndCorrelations

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = MAPFromFewMomentsAndCorrelations(moms, corr1, r)`
        * - Mathematica:
          - :code:`{D0, D1} = MAPFromFewMomentsAndCorrelations[moms, corr1, r]`
        * - Python/Numpy:
          - :code:`D0, D1 = MAPFromFewMomentsAndCorrelations(moms, corr1, r)`

    Creates a Markovian arrival process that has the given
    2 or 3 marginal moments and lag-1 autocorrelation.
    The decay of the autocorrelation function can be optionally
    adjusted as well.
    The lag-k autocorrelation function :math:`\rho_k` of the 
    resulting MAP is :math:`\rho_k=r(corr_1/r)^k`.

    Parameters
    ----------
    moms : vector of doubles, length 2 or 3
        The list of marginal moments to match. 
    corr1 : double
        The lag-1 autocorrelation coefficient to match.
    r : double, optional
        The decay of the autocorrelation function.
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    
    Notes
    -----
    With 2 marginal moments, or with 3 marginal moments and 
    positive autocorrelation the procedure always returns a 
    valid Markovian representation.

    References
    ----------
    .. [1] G Horvath, "Matching marginal moments and lag 
           autocorrelations with MAPs," ValueTools 2013, 
           Torino, Italy (2013).

    Examples
    ========
    For Matlab:

    >>> moms = [1.1, 6.05];
    >>> corr1 = -0.17;
    >>> [D0, D1] = MAPFromFewMomentsAndCorrelations(moms, corr1);
    >>> disp(D0);
         -0.28494      0.28494            0            0
                0      -18.134            0            0
                0            0     -0.28494      0.28494
                0            0            0     -0.95707
    >>> disp(D1);
                0            0            0            0
         0.022795       1.4279       4.9669       11.716
                0            0            0            0
         0.013835      0.86667     0.022795     0.053771
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 2);
    >>> disp(rmoms);
              1.1         6.05
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1);
    >>> disp(rcorr1);
            -0.17
    >>> moms = [1.2, 4.32, 20.];
    >>> corr1 = -0.4;
    >>> [D0, D1] = MAPFromFewMomentsAndCorrelations(moms, corr1);
    >>> disp(D0);
         -0.33604      0.33604            0            0            0            0
                0      -36.027            0            0            0            0
                0            0      -1.2829       1.2829            0            0
                0            0            0      -1.2829       1.2829            0
                0            0            0            0      -1.2829       1.2829
                0            0            0            0            0      -1.3747
    >>> disp(D1);
                0            0            0            0            0            0
         0.024741       1.7766       23.475            0            0       10.751
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
         0.017937        1.288     0.047145            0            0     0.021591
    >>> BuToolsCheckPrecision = 10.^-12;
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(rmoms);
              1.2         4.32           20
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1);
    >>> disp(rcorr1);
             -0.4
    >>> moms = [1.2, 4.32, 20.];
    >>> corr1 = 0.4;
    >>> [D0, D1] = MAPFromFewMomentsAndCorrelations(moms, corr1);
    >>> disp(D0);
      Columns 1 through 6
         -0.54184      0.54184            0            0            0            0
                0      -116.34            0            0            0            0
                0            0     -0.24417      0.24417            0            0
                0            0            0       -2.014        2.014            0
                0            0            0            0       -2.014        2.014
                0            0            0            0            0       -2.014
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
      Columns 7 through 9
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
            2.014            0            0
           -2.014        2.014            0
                0       -2.014        2.014
                0            0       -2.014
    >>> disp(D1);
      Columns 1 through 6
                0            0            0            0            0            0
           16.128       88.575     0.099774       11.534            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
                0            0            0            0            0            0
         0.062048      0.34076     0.013818       1.5974            0            0
      Columns 7 through 9
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
                0            0            0
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(rmoms);
              1.2         4.32           20
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1);
    >>> disp(rcorr1);
              0.4

    For Mathematica:

    
    For Python/Numpy:

    >>> moms = [1.1, 6.05]
    >>> corr1 = -0.17
    >>> D0, D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
    >>> print(D0)
    [[ -0.28494   0.28494   0.        0.     ]
     [ -0.      -18.13384   0.        0.     ]
     [  0.        0.       -0.28494   0.28494]
     [  0.        0.       -0.       -0.95707]]
    >>> print(D1)
    [[  0.        0.        0.        0.     ]
     [  0.0228    1.42791   4.9669   11.71623]
     [  0.        0.        0.        0.     ]
     [  0.01384   0.86667   0.0228    0.05377]]
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 2)
    >>> print(rmoms)
    [1.0999999999999996, 6.0499999999999963]
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1)[0]
    >>> print(rcorr1)
    -0.17
    >>> moms = [1.2, 4.32, 20.]
    >>> corr1 = -0.4
    >>> D0, D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
    >>> print(D0)
    [[ -0.33604   0.33604   0.        0.        0.        0.     ]
     [  0.      -36.02667   0.        0.        0.        0.     ]
     [  0.        0.       -1.28286   1.28286   0.        0.     ]
     [  0.        0.        0.       -1.28286   1.28286   0.     ]
     [  0.        0.        0.        0.       -1.28286   1.28286]
     [  0.        0.        0.        0.        0.       -1.37472]]
    >>> print(D1)
    [[  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  2.47407e-02   1.77659e+00   2.34746e+01   0.00000e+00   0.00000e+00   1.07508e+01]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  1.79373e-02   1.28805e+00   4.71448e-02   0.00000e+00   0.00000e+00   2.15912e-02]]
    >>> butools.checkPrecision = 10.**-12
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(rmoms)
    [1.2000000000000002, 4.3200000000000021, 19.999999999999915]
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1)[0]
    >>> print(rcorr1)
    -0.4
    >>> moms = [1.2, 4.32, 20.]
    >>> corr1 = 0.4
    >>> D0, D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
    >>> print(D0)
    [[  -0.54184    0.54184    0.         0.         0.         0.         0.         0.         0.     ]
     [   0.      -116.33647    0.         0.         0.         0.         0.         0.         0.     ]
     [   0.         0.        -0.24417    0.24417    0.         0.         0.         0.         0.     ]
     [   0.         0.         0.        -2.01403    2.01403    0.         0.         0.         0.     ]
     [   0.         0.         0.         0.        -2.01403    2.01403    0.         0.         0.     ]
     [   0.         0.         0.         0.         0.        -2.01403    2.01403    0.         0.     ]
     [   0.         0.         0.         0.         0.         0.        -2.01403    2.01403    0.     ]
     [   0.         0.         0.         0.         0.         0.         0.        -2.01403    2.01403]
     [   0.         0.         0.         0.         0.         0.         0.         0.        -2.01403]]
    >>> print(D1)
    [[  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  1.61283e+01   8.85745e+01   9.97739e-02   1.15339e+01   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  6.20477e-02   3.40758e-01   1.38184e-02   1.59740e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]]
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(rmoms)
    [1.2, 4.3199999999999985, 19.999999999999851]
    >>> rcorr1 = LagCorrelationsFromMAP(D0, D1, 1)[0]
    >>> print(rcorr1)
    0.4

