butools.ph.MEFromMoments
========================

.. currentmodule:: butools.ph

.. np:function:: MEFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MEFromMoments(moms)`
        * - Mathematica:
          - :code:`{alpha, A} = MEFromMoments[moms]`
        * - Python/Numpy:
          - :code:`alpha, A = MEFromMoments(moms)`

    Creates a matrix-exponential distribution that has the
    same moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(2*M-1)
        The list of moments. The order of the resulting 
        matrix-exponential distribution is 
        determined based on the number of moments given. To 
        obtain a matrix exponential distribution of order M,
        2*M-1 moments are required.

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.

    References
    ----------
    .. [1] A. van de Liefvoort. The moment problem for 
           continuous distributions. Technical report, 
           University of Missouri, WP-CM-1990-02, Kansas City,
           1990.

    Examples
    ========
    For Matlab:

    >>> a = [0.1, 0.9, 0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> moms = MomentsFromPH(a,A,5);
    >>> disp(moms);
          0.20939      0.10449     0.089092      0.11027      0.17953
    >>> [a,A] = MEFromMoments(moms);
    >>> disp(a);
          0.33333      0.33333      0.33333
    >>> disp(A);
           -8.085       10.177      -9.9999
          -1.2584      -5.1438       1.7873
          -1.9255       1.9599      -4.9712
    >>> memoms = MomentsFromME(a,A,5);
    >>> disp(memoms);
          0.20939      0.10449     0.089092      0.11027      0.17953

    For Mathematica:

    >>> a = {0.1, 0.9, 0};
    >>> A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}};
    >>> moms = MomentsFromPH[a,A,5];
    >>> Print[moms];
    {0.20938722294654497, 0.10448912014333092, 0.08909150039190732, 0.11026674096545433, 0.179530273247209}
    >>> {a,A} = MEFromMoments[moms];
    >>> Print[a];
    {1/3, 1/3, 1/3}
    >>> Print[A];
    {{-8.084991901863543, 10.176612950254244, -9.999939852873878},
     {-1.2584436954049303, -5.143808678565172, 1.7873056877319793},
     {-1.925450036434008, 1.959914946726235, -4.971199419570917}}
    >>> memoms = MomentsFromME[a,A,5];
    >>> Print[memoms];
    {0.209387222946545, 0.10448912014333095, 0.08909150039190736, 0.11026674096545445, 0.17953027324720922}

    For Python/Numpy:

    >>> a = ml.matrix([[0.1, 0.9, 0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> moms = MomentsFromPH(a,A,5)
    >>> print(moms)
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288, 0.11026674096545433, 0.17953027324720897]
    >>> a,A = MEFromMoments(moms)
    >>> print(a)
    [[ 0.33333  0.33333  0.33333]]
    >>> print(A)
    [[ -8.08499  10.17661  -9.99994]
     [ -1.25844  -5.14381   1.78731]
     [ -1.92545   1.95991  -4.9712 ]]
    >>> memoms = MomentsFromME(a,A,5)
    >>> print(memoms)
    [0.20938722294654494, 0.10448912014333089, 0.089091500391907275, 0.11026674096545437, 0.1795302732472091]

