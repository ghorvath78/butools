butools.fitting.RelativeEntropy
===============================

.. currentmodule:: butools.fitting

.. np:function:: RelativeEntropy

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`re = RelativeEntropy(p1, p2)`
        * - Mathematica:
          - :code:`re = RelativeEntropy[p1, p2]`
        * - Python/Numpy:
          - :code:`re = RelativeEntropy(p1, p2)`

    Returns the relative entropy (aka Kullback–Leibler
    divergence) of two vectors.

    Parameters
    ----------
    p1 : vector, length M
        The first vector
    p2 : vector, length M
        The second vector

    Returns
    -------
    re : double
        The relative entropy calculated as
        :math:`re=\sum_{i=1}^M p1_i |\log(p1_i/p2_i)|`

    Examples
    ========
    For Matlab:

    >>> tr = dlmread('/home/gabor/github/butools/test/data/bctrace.iat');
    >>> trAcf = LagCorrelationsFromTrace(tr, 10);
    >>> disp(trAcf);
      Columns 1 through 6
          0.20005      0.18927      0.13895      0.14213      0.11713      0.12368
      Columns 7 through 10
          0.11212      0.10051      0.10019     0.098797
    >>> [D0, D1] = MAPFromFewMomentsAndCorrelations(MarginalMomentsFromTrace(tr, 3), trAcf(1));
    >>> mapAcf = LagCorrelationsFromMAP(D0, D1, 10);
    >>> disp(mapAcf);
      Columns 1 through 6
          0.20005      0.12003     0.072023     0.043216     0.025931     0.015559
      Columns 7 through 10
        0.0093357    0.0056017    0.0033611    0.0020168
    >>> reAcf = RelativeEntropy(mapAcf, trAcf);
    >>> disp(reAcf);
          0.28344

    For Mathematica:

    >>> tr = Flatten[Import["/home/gabor/github/butools/test/data/bctrace.iat","CSV"]];
    >>> trAcf = LagCorrelationsFromTrace[tr, 10];
    >>> Print[trAcf];
    {0.2000486547928462, 0.18927461196895048, 0.1389510350355339, 0.14213385937998096, 0.11712708793961699, 0.12367812078281121, 0.11212100297743381, 0.10051058879698098, 0.10019060948165948, 0.0987971115012499}
    >>> {D0, D1} = MAPFromFewMomentsAndCorrelations[MarginalMomentsFromTrace[tr, 3], trAcf[[1]]];
    >>> mapAcf = LagCorrelationsFromMAP[D0, D1, 10];
    >>> Print[mapAcf];
    {0.20004865479284564, 0.12003405953863643, 0.07202335583933256, 0.04321576564432875, 0.025930510713659746, 0.015558937250009402, 0.009335740858440003, 0.005601671629332731, 0.0033611392516858904, 0.0020167653187785875}
    >>> reAcf = RelativeEntropy[mapAcf, trAcf];
    >>> Print[reAcf];
    0.2834377149084656

    For Python/Numpy:

    >>> tr = np.loadtxt("/home/gabor/github/butools/test/data/bctrace.iat")
    >>> trAcf = LagCorrelationsFromTrace(tr, 10)
    >>> print(trAcf)
    [0.20004885484210852, 0.18927480124417365, 0.13895117398714094, 0.14213400151438738, 0.11712720506721902, 0.12367824446146115, 0.1121211150989561, 0.10051068930807179, 0.10019070967277495, 0.098797210298893726]
    >>> D0, D1 = MAPFromFewMomentsAndCorrelations(MarginalMomentsFromTrace(tr, 3), trAcf[0])
    >>> mapAcf = LagCorrelationsFromMAP(D0, D1, 10)
    >>> print(mapAcf)
    [ 0.20005  0.12003  0.07202  0.04322  0.02593  0.01556  0.00934  0.0056   0.00336  0.00202]
    >>> reAcf = RelativeEntropy(mapAcf, trAcf)
    >>> print(reAcf)
    0.283438051725

