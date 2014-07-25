butools.dmap.ImageFromDMMAP
===========================

.. currentmodule:: butools.dmap

.. np:function:: ImageFromDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`ImageFromDMMAP(D, outFileName, prec)`
        * - Mathematica:
          - :code:`ImageFromDMMAP[D, outFileName, prec]`
        * - Python/Numpy:
          - :code:`ImageFromDMMAP(D, outFileName, prec)`

    Depicts the given discrete marked Markovian arrival
    process, and either displays it or saves it to file.

    Parameters
    ----------
    D : list of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
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

    >>> D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
    >>> D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
    >>> D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
    >>> D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]
    >>> D = {D0,D1,D2,D3};
    >>> ImageFromDMMAP(D,'figure.pdf');
    >>> ImageFromDMMAP(D);

