butools.map.ImageFromMMAP
=========================

.. currentmodule:: butools.map

.. np:function:: ImageFromMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`ImageFromMMAP(D, outFileName, prec)`
        * - Mathematica:
          - :code:`ImageFromMMAP[D, outFileName, prec]`
        * - Python/Numpy:
          - :code:`ImageFromMMAP(D, outFileName, prec)`

    Depicts the given marked Markovian arrival process, and
    either displays it or saves it to file.

    Parameters
    ----------
    D : list of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
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

    >>> D0=[-1.78 0.29; 0.07 -0.92]
    >>> D1=[0.15 0.49; 0.23 0.36]
    >>> D2=[0.11 0.2; 0.01 0]
    >>> D3=[0.14 0.4; 0.11 0.14]
    >>> D = {D0,D1,D2,D3};
    >>> ImageFromMMAP(D,'figure.pdf');
    >>> ImageFromMMAP(D);

