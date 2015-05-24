butools.map.ImageFromMAP
========================

.. currentmodule:: butools.map

.. np:function:: ImageFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`ImageFromMAP(D0, D1, outFileName, prec)`
        * - Mathematica:
          - :code:`ImageFromMAP[D0, D1, outFileName, prec]`
        * - Python/Numpy:
          - :code:`ImageFromMAP(D0, D1, outFileName, prec)`

    Depicts the given Markovian arrival process, and either
    displays it or saves it to file.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    >>> D0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9];
    >>> D1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1];
    >>> ImageFromMAP(D0,D1,'figure.pdf');
    >>> ImageFromMAP(D0,D1);

