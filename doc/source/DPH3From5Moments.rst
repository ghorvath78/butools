butools.dph.DPH3From5Moments
============================

.. currentmodule:: butools.dph

.. np:function:: DPH3From5Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = DPH3From5Moments(moms, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = DPH3From5Moments[moms, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = DPH3From5Moments(moms, prec)`

    Returns an order-3 discrete phase-type distribution 
    which has the same 5 moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(5)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,3)
      The initial probability vector of the DPH(3)
    A : matrix, shape (3,3)
      Transition probability matrix of the DPH(3)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a DPH(3).

    This procedure first calls 'MGFromMoments', then transforms
    it to DPH(3) by 'CanonicalFromDPH3'.

    Examples
    ========
    For Matlab:

    >>> a = [0.7,0.1,0.2];
    >>> A = [0.2, 0.51, 0.1; 0.58, 0.41, 0; 0.1, 0.4, 0.3];
    >>> moms = MomentsFromDPH(a, A);
    >>> disp(moms);
           9.3096        175.1       4968.7   1.8805e+05   8.8966e+06
    >>> [b, B] = DPH3From5Moments(moms);
    >>> disp(b);
          0.73989     0.076837      0.18327
    >>> disp(B);
          0.89971      0.10029            0
                0     0.010293      0.98971
                0     0.050581            0
    >>> phmoms = MomentsFromMG(b, B, 5);
    >>> disp(phmoms);
           9.3096        175.1       4968.7   1.8805e+05   8.8966e+06

    For Mathematica:

    >>> a = {0.7,0.1,0.2};
    >>> A = {{0.2, 0.51, 0.1},{0.58, 0.41, 0},{0.1, 0.4, 0.3}};
    >>> moms = MomentsFromDPH[a, A];
    >>> Print[moms];
    {9.3096349745331, 175.10327171027384, 4968.663522150066, 188050.43861214988, 8.896632715174045*^6}
    >>> {b, B} = DPH3From5Moments[moms];
    "Ordered eigenvalues:"{0.8997069306212445, 0.22894830561223567, -0.21865523626858424}
    >>> Print[b];
    {0.7398925149830309, 0.07683743218204545, 0.18327005283492365}
    >>> Print[B];
    {{0.8997069306212445, 0.10029306937875548, 0},
     {0, 0.010293069343651429, 0.9897069306563486},
     {0, 0.05058138354526465, 0}}
    >>> phmoms = MomentsFromMG[b, B, 5];
    >>> Print[phmoms];
    {9.309634974533111, 175.10327171027419, 4968.66352215008, 188050.43861215067, 8.89663271517409*^6}

    For Python/Numpy:

    >>> a = ml.matrix([[0.7,0.1,0.2]])
    >>> A = ml.matrix([[0.2, 0.51, 0.1],[0.58, 0.41, 0],[0.1, 0.4, 0.3]])
    >>> moms = MomentsFromDPH(a, A)
    >>> print(moms)
    [9.3096349745331022, 175.10327171027393, 4968.6635221500701, 188050.43861215009, 8896632.7151740566]
    >>> b, B = DPH3From5Moments(moms)
    >>> print(b)
    [[ 0.73989  0.07684  0.18327]]
    >>> print(B)
    [[ 0.89971  0.10029  0.     ]
     [ 0.       0.01029  0.98971]
     [ 0.       0.05058  0.     ]]
    >>> phmoms = MomentsFromMG(b, B, 5)
    >>> print(phmoms)
    [9.3096349745330773, 175.10327171027262, 4968.6635221500082, 188050.4386121468, 8896632.7151738554]

