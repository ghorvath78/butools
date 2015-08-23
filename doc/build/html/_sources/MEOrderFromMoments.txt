butools.ph.MEOrderFromMoments
=============================

.. currentmodule:: butools.ph

.. np:function:: MEOrderFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`order = MEOrderFromMoments(moms, prec)`
        * - Mathematica:
          - :code:`order = MEOrderFromMoments[moms, prec]`
        * - Python/Numpy:
          - :code:`order = MEOrderFromMoments(moms, prec)`

    Returns the order of ME distribution that can realize
    the given moments.
    
    Parameters
    ----------
    moms : list of doubles
        The list of moments
    prec : double, optional
        Precision used to detect if the determinant of the
        Hankel matrix is zero. The default value is 1e-12.
    
    Returns
    -------
    order : int
        The order of ME distribution that can realize the 
        given moments

    References
    ----------
    .. [1]  L. Bodrog, A. Horvath, M. Telek, "Moment 
            characterization of matrix exponential and Markovian
            arrival processes," Annals of Operations Research, 
            vol. 160, pp. 51-68, 2008.

    Examples
    ========
    For Matlab:

    >>> a = [0.1,0.9,0];
    >>> A = [-6.2, 2., 0.; 2., -9., 1.; 1., 0., -3.];
    >>> moms = MomentsFromME(a, A);
    >>> disp(moms);
          0.20939      0.10449     0.089092      0.11027      0.17953
    >>> mo = MEOrderFromMoments(moms);
    >>> disp(mo);
         3
    >>> b = [0.2,0.3,0.5];
    >>> B = [-1., 0., 0.; 0., -3., 2.; 0., -2., -3.];
    >>> [a, A] = MonocyclicPHFromME(b, B);
    >>> moms = MomentsFromME(a, A);
    >>> disp(moms);
      Columns 1 through 6
          0.35385      0.41893       1.1552       4.6998       23.838       143.78
      Columns 7 through 12
           1007.8       8064.3        72578   7.2577e+05   7.9834e+06     9.58e+07
      Columns 13 through 18
       1.2454e+09   1.7436e+10   2.6153e+11   4.1846e+12   7.1137e+13   1.2805e+15
      Columns 19 through 24
       2.4329e+16   4.8658e+17   1.0218e+19    2.248e+20   5.1704e+21   1.2409e+23
      Columns 25 through 30
       3.1022e+24   8.0658e+25   2.1778e+27   6.0978e+28   1.7684e+30   5.3051e+31
      Columns 31 through 36
       1.6446e+33   5.2626e+34   1.7367e+36   5.9047e+37   2.0666e+39   7.4399e+40
      Columns 37 through 42
       2.7528e+42    1.046e+44   4.0796e+45   1.6318e+47   6.6905e+48     2.81e+50
      Columns 43 through 48
       1.2083e+52   5.3165e+53   2.3924e+55   1.1005e+57   5.1725e+58   2.4828e+60
      Columns 49 through 53
       1.2166e+62   6.0828e+63   3.1022e+65   1.6132e+67   8.5498e+68
    >>> mo = MEOrderFromMoments(moms);
    >>> disp(mo);
         3

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[0.1,0.9,0]])
    >>> A = ml.matrix([[-6.2, 2., 0.],[2., -9., 1.],[1., 0., -3.]])
    >>> moms = MomentsFromME(a, A)
    >>> print(moms)
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288, 0.11026674096545433, 0.17953027324720897]
    >>> mo = MEOrderFromMoments(moms)
    >>> print(mo)
    3
    >>> b = ml.matrix([[0.2,0.3,0.5]])
    >>> B = ml.matrix([[-1., 0., 0.],[0., -3., 2.],[0., -2., -3.]])
    >>> a, A = MonocyclicPHFromME(b, B)
    >>> moms = MomentsFromME(a, A)
    >>> print(moms)
    [0.35384615384615531, 0.41893491124260573, 1.155211652253076, 4.699835439935578, 23.83775616561579, 143.78185836646887, 1007.8194071104439, 8064.272882521419, 72578.133718784462, 725767.95874615829, 7983382.3513983944, 95800362.980475202, 1245404149.6660547, 17435657571.499924, 261534869138.85153, 4184557956997.604, 71137485550890.109, 1280474741063026.2, 24329020082838864.0, 4.8658040164763066e+17, 1.021818843442624e+19, 2.248001455559353e+20, 5.17040334777797e+21, 1.2408968034663769e+23, 3.1022420086659455e+24, 8.065829222531619e+25, 2.1777738900835537e+27, 6.0977668922339555e+28, 1.7683523987478442e+30, 5.305057196243521e+31, 1.6445677308354878e+33, 5.262616738673549e+34, 1.7366635237622677e+36, 5.904655980791697e+37, 2.066629593277089e+39, 7.439866535797508e+40, 2.752750618245071e+42, 1.0460452349331244e+44, 4.0795764162391765e+45, 1.6318305664956671e+47, 6.690505322632222e+48, 2.810012235505527e+50, 1.2083052612673737e+52, 5.316543149576433e+53, 2.392444417309389e+55, 1.1005244319623167e+57, 5.172464830222877e+58, 2.4827831185069758e+60, 1.2165637280684155e+62, 6.082818640342062e+63, 3.1022375065744455e+65, 1.613163503418708e+67, 8.549766568119135e+68]
    >>> mo = MEOrderFromMoments(moms)
    >>> print(mo)
    3

