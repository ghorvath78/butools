butools.dmap.ImageFromDMAP
==========================

.. currentmodule:: butools.dmap

.. np:function:: ImageFromDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`ImageFromDMAP(D0, D1, outFileName, prec)`
        * - Mathematica:
          - :code:`ImageFromDMAP[D0, D1, outFileName, prec]`
        * - Python/Numpy:
          - :code:`ImageFromDMAP(D0, D1, outFileName, prec)`

    Depicts the given discrete Markovian arrival process,
    and either displays it or saves it to file.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete MAP.
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete MAP.
    outFileName : string, optional
        If it is not provided, or equals to 'display', the
        image is displayed on the screen, otherwise it is 
        written to the file. The file format is deduced 
        from the file name.
    prec : double, optional
        Transition probabilities less then prec are 
        considered to be zero and are left out from the 
        image. The default value is 1e-13.

    Notes
    -----
    The 'graphviz' software must be installed and available
    in the path to use this feature.

    Examples
    --------
    For Matlab:

    >>> D0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.24];
    >>> D1=[0 0.88 0.1; 0.42 0.07 0.14; 0.13 0.15 0.15];
    >>> ImageFromDMAP(D0,D1,'figure.pdf');
    >>> ImageFromDMAP(D0,D1);

