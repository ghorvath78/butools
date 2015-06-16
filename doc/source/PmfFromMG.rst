butools.dph.PmfFromMG
=====================

.. currentmodule:: butools.dph

.. np:function:: PmfFromMG

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pmf = PmfFromMG(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`pmf = PmfFromMG[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`pmf = PmfFromMG(alpha, A, x, prec)`

    Returns the probability mass function of a matrix-
    geometric distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric
        distribution. The sum of the entries of pi0 is less
        or equal to 1.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric
        distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input MG 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    pmf : column vector of doubles
        The probabilities that the matrix-geometrically 
        distributed random variable takes the corresponding "x"
        values
        
    Examples
    ========
    For Matlab:

    >>> a = [-0.6, 0.3, 1.3];
    >>> A = [0.25, 0.2, -0.15; 0.3, 0.1, 0.25; 0, 0.2, 0.47];
    >>> x = (0:1:100);
    >>> pmf = PmfFromMG(a,A,x);
    >>> plot(x,pmf);

    For Mathematica:

    >>> a = {-0.6, 0.3, 1.3};
    >>> A = {{0.25, 0.2, -0.15},{0.3, 0.1, 0.25},{0, 0.2, 0.47}};
    >>> x = Range[0,100,1];
    >>> pmf = PmfFromMG[a,A,x];
    >>> ListLinePlot[Transpose[{x, pmf}]]

    For Python/Numpy:

    >>> a = ml.matrix([[-0.6, 0.3, 1.3]])
    >>> A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    >>> x = np.arange(0,101.0,1)
    >>> pmf = PmfFromMG(a,A,x)
    >>> plt.plot(x,pmf)

