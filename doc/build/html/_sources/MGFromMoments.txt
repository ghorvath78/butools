butools.dph.MGFromMoments
=========================

.. currentmodule:: butools.dph

.. np:function:: MGFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MGFromMoments(moms)`
        * - Mathematica:
          - :code:`{alpha, A} = MGFromMoments[moms]`
        * - Python/Numpy:
          - :code:`alpha, A = MGFromMoments(moms)`

    Creates a matrix-geometric distribution that has the
    same moments as given.

    Parameters
    ----------
    moms : vector of doubles
        The list of moments. The order of the resulting 
        matrix-geometric distribution is 
        determined based on the number of moments given. To 
        obtain a matrix-geometric distribution of order M,
        2*M-1 moments are required.

    Returns
    -------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
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

    >>> moms = [4.08, 20.41, 130.45, 1054.41, 10463.73];
    >>> [a,A] = MGFromMoments(moms);
    >>> disp(a);
          0.33333      0.33333      0.33333
    >>> disp(A);
          0.15523       1.7289      0.10482
        -0.013774       0.6823    -0.023472
        -0.013847     -0.16787      0.82688
    >>> memoms = MomentsFromMG(a,A,5);
    >>> disp(memoms);
             4.08        20.41       130.45       1054.4        10464

    For Mathematica:

    >>> moms = {4.08, 20.41, 130.45, 1054.41, 10463.73};
    >>> {a,A} = MGFromMoments[moms];
    >>> Print[a];
    {1/3, 1/3, 1/3}
    >>> Print[A];
    {{0.15522721633282086, 1.7288735256877237, 0.10482133882430097},
     {-0.013773788451490479, 0.6823009288291466, -0.02347241196722473},
     {-0.013846712675345957, -0.1678656131152182, 0.8268849851301606}}
    >>> memoms = MomentsFromMG[a,A,5];
    >>> Print[memoms];
    {4.080000000000002, 20.41000000000002, 130.4500000000002, 1054.4100000000037, 10463.730000000038}

    For Python/Numpy:

    >>> moms = [4.08, 20.41, 130.45, 1054.41, 10463.73]
    >>> a,A = MGFromMoments(moms)
    >>> print(a)
    [[ 0.33333  0.33333  0.33333]]
    >>> print(A)
    [[ 0.15523  1.72887  0.10482]
     [-0.01377  0.6823  -0.02347]
     [-0.01385 -0.16787  0.82688]]
    >>> memoms = MomentsFromMG(a,A,5)
    >>> print(memoms)
    [4.080000000000001, 20.410000000000029, 130.45000000000047, 1054.4100000000135, 10463.730000000427]

