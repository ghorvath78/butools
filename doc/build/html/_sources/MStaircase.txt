butools.reptrans.MStaircase
===========================

.. currentmodule:: butools.reptrans

.. np:function:: MStaircase

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[B, n] = MStaircase(X, Z, precision)`
        * - Mathematica:
          - :code:`{B, n} = MStaircase[X, Z, precision]`
        * - Python/Numpy:
          - :code:`B, n = MStaircase(X, Z, precision)`

    Computes a smaller representation using the staircase 
    algorithm.
    
    Notes
    -----
    This function should not be called directly.
    It is used by 'MinimalRepFromME' and 'MinimalRepFromRAP'.
    
    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation 
            of rational arrival processes." Madrid Conference
            on Qeueuing theory (MCQT), June 2010.

