butools.dph.DPH3From5Moments
============================

.. currentmodule:: butools.dph

.. np:function:: DPH3From5Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = DPH3From5Moments(moms, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = DPH3From5Moments[moms, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = DPH3From5Moments(moms, prec)`

    Returns an order-3 discrete phase-type distribution 
    which has the same 5 moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(5)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,3)
      The initial probability vector of the DPH(3)
    A : matrix, shape (3,3)
      Transition probability matrix of the DPH(3)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a DPH(3).

    This procedure first calls 'MGFromMoments', then transforms
    it to DPH(3) by 'CanonicalFromDPH3'.

    Examples
    --------
    For Matlab:
    
    >>> moms = [9.3096, 175.1, 4968.7, 1.8805e+05, 8.8966e+06];
    >>> [a,A]=DPH3From5Moments(moms)
    >>> a
      0.73989     0.076837      0.18327
    >>> A
      0.89971      0.10029            0
            0     0.010293      0.98971
            0     0.050581            0
    >>> MomentsFromDPH(a,A)
       9.3096        175.1       4968.7   1.8805e+05   8.8966e+06

    For Python/Numpy:
    
    >>> moms = [9.3096349745331022, 175.10327171027393, 4968.6635221500701, 188050.43861215009, 8896632.7151740566]
    >>> a,A=DPH3From5Moments(moms)
    >>> print(a)
    [[ 0.73989251  0.07683743  0.18327005]]
    >>> print(A)
    [[ 0.89970693  0.10029307  0.        ]
     [ 0.          0.01029307  0.98970693]
     [ 0.          0.05058138  0.        ]]
    >>> print(MomentsFromDPH(a,A))
    [9.3096349745330773, 175.10327171027262, 4968.6635221500082, 188050.4386121468, 8896632.7151738554]
    
