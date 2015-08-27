butools.ph.PdfFromPH
====================

.. currentmodule:: butools.ph

.. np:function:: PdfFromPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pdf = PdfFromPH(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`pdf = PdfFromPH[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`pdf = PdfFromPH(alpha, A, x, prec)`

    Returns the probability density function of a continuous
    phase-type distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    x : vector of doubles
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input ME 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    pdf : column vector of doubles
        The values of the density function at the 
        corresponding "x" values

    Examples
    ========
    For Matlab:

    >>> a = [0.1,0.9,0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> x = (0.0:0.002:3.0);
    >>> pdf = PdfFromPH(a, A, x);
    >>> plot(x, pdf)

    For Mathematica:

    >>> a = {0.1,0.9,0};
    >>> A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}};
    >>> x = Range[0.0,3.0,0.002];
    >>> pdf = PdfFromPH[a, A, x];
    >>> ListLinePlot[{Transpose[{x, pdf}]}]

    For Python/Numpy:

    >>> a = ml.matrix([[0.1,0.9,0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> x = np.arange(0.0,3.002,0.002)
    >>> pdf = PdfFromPH(a, A, x)
    >>> plt.plot(x, pdf)

