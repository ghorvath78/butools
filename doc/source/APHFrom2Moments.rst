butools.ph.APHFrom2Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: APHFrom2Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = APHFrom2Moments(moms, maxSize)`
        * - Mathematica:
          - :code:`{alpha, A} = APHFrom2Moments[moms, maxSize]`
        * - Python/Numpy:
          - :code:`alpha, A = APHFrom2Moments(moms, maxSize)`

    Returns an acyclic PH which has the same 2 moments as
    given. If detects the order and the structure 
    automatically to match the given moments.

    Parameters
    ----------
    moms : vector of doubles, length(2)
      The moments to match
    maxSize : int, optional
      The maximal size of the resulting APH. The default value
      is 100.

    Returns
    -------
    alpha : vector, shape (1,M)
      The initial probability vector of the APH
    A : matrix, shape (M,M)
      Transient generator matrix of the APH
    
    Raises an error if the moments are not feasible with an
    APH of size "maxSize".
    
    Examples
    --------

