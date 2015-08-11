butools.reptrans.FindMarkovianRepresentation
============================================

.. currentmodule:: butools.reptrans

.. np:function:: FindMarkovianRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`mrep = FindMarkovianRepresentation(rep, @transfun, @evalfunc, precision)`
        * - Mathematica:
          - :code:`mrep = FindMarkovianRepresentation[rep, transfun, evalfunc, precision]`
        * - Python/Numpy:
          - :code:`mrep = FindMarkovianRepresentation(rep, transfun, evalfunc, precision)`
    
    Obtains a Markovian representation from a non-Markovian 
    one while keeping the size the same, by applying a series 
    of elementary transformations.
    
    Parameters
    ----------
    rep : tuple of matrices
        The initial non-Markovian representation
        (initial vector and generator of a PH, matrices of a
        MAP, or a MMAP, etc.)
    transfun : callable
        A function that transforms the representation using 
        the given similarity transformation matrix
    evalfunc : callable
        A function that returns how far the representation is
        from the Markovian one
    precision : double
        A representation is considered to be a Markovian one
        if it is closer than the precision. The default value
        is 1e-7
        
    Returns
    -------
    mrep : tuple of matrices
        The Markovian representation, if found. If not found,
        the closest one is returned.

    Notes
    -----
    This function should not be called directly.
    It is used by 'PHFromME', 'MAPFromRAP', etc. functions.

    References
    ----------
    .. [1]  G Horváth, M Telek, "A minimal representation of 
            Markov arrival processes and a moments matching 
            method," Performance Evaluation 64:(9-12) 
            pp. 1153-1168. (2007)

