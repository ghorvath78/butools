butools.ph.ImageFromPH
======================

.. currentmodule:: butools.ph

.. np:function:: ImageFromPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`ImageFromPH(alpha, A, outFileName, prec)`
        * - Mathematica:
          - :code:`ImageFromPH[alpha, A, outFileName, prec]`
        * - Python/Numpy:
          - :code:`ImageFromPH(alpha, A, outFileName, prec)`

    Depicts the given phase-type distribution,
    and either displays it or saves it to file.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    outFileName : string, optional
        If it is not provided, or equals to 'display', the
        image is displayed on the screen, otherwise it is 
        written to the file. The file format is deduced 
        from the file name.
    prec : double, optional
        Transition rates less then prec are considered to
        be zero and are left out from the image. The 
        default value is 1e-13.

    Notes
    -----
    The 'graphviz' software must be installed and available
    in the path to use this feature.

    Examples
    --------
    For Matlab:

    >>> [alpha,A]=RandomPH(3,8,7);
    >>> ImageFromPH(alpha,A,'figure.pdf');
    >>> ImageFromPH(alpha,A);

