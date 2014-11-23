butools.ph.PdfFromME
====================

.. currentmodule:: butools.ph

.. np:function:: PdfFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pdf = PdfFromME(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`pdf = PdfFromME[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`pdf = PdfFromME(alpha, A, x, prec)`

    Returns the probability density function of a matrix-
    exponential distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
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
    --------    
    For Matlab:
    
    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> x = (0:0.01:2)';
    >>> pdf = PdfFromME(a, A, x);
    >>> plot(x, pdf)

    For Python/Numpy:
    
    >>> a = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    >>> x = np.linspace(0,2,201)
    >>> pdf = PdfFromME(a,A,x)
    >>> plt.plot(x,pdf)

