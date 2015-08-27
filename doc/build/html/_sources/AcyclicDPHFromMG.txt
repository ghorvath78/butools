butools.dph.AcyclicDPHFromMG
============================

.. currentmodule:: butools.dph

.. np:function:: AcyclicDPHFromMG

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = AcyclicDPHFromMG(alpha, A, precision)`
        * - Mathematica:
          - :code:`{beta, B} = AcyclicDPHFromMG[alpha, A, precision]`
        * - Python/Numpy:
          - :code:`beta, B = AcyclicDPHFromMG(alpha, A, precision)`

    Transforms a matrix-geometric representation to an acyclic
    DPH representation of the same size, if possible.
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the acyclic discrete
        phase-type representation
    B : matrix, shape (M,M)
        Transition probability matrix of the acyclic discrete
        phase-type representation
    
    Notes
    -----
    Contrary to 'AcyclicPHFromME' of the 'ph' package, this 
    procedure is not able to extend the size in order to obtain
    a Markovian initial vector.

    Raises an error if A has complex eigenvalues. In this case
    the transformation to an acyclic representation is not 
    possible

    Examples
    ========
    For Matlab:

    >>> a = [0,0,1.0];
    >>> A = [0.22, 0, 0; 0.3, 0.1, 0.55; 0.26, 0, 0.73];
    >>> [b, B] = AcyclicDPHFromMG(a, A);
    >>> disp(b);
          0.69103      0.29786     0.011111
    >>> disp(B);
             0.73         0.27            0
                0         0.22         0.78
                0            0          0.1
    >>> ma = MomentsFromMG(a, A, 5);
    >>> disp(ma);
           4.9383       34.807       339.49       4335.8        68954
    >>> mb = MomentsFromMG(b, B, 5);
    >>> disp(mb);
           4.9383       34.807       339.49       4335.8        68954
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
        2.999e-16

    For Mathematica:

    >>> a = {0,0,1.0};
    >>> A = {{0.22, 0, 0},{0.3, 0.1, 0.55},{0.26, 0, 0.73}};
    >>> {b, B} = AcyclicDPHFromMG[a, A];
    >>> Print[b];
    {0.691025641025641, 0.2978632478632479, 0.011111111111111127}
    >>> Print[B];
    {{0.73, 0.27, 0.},
     {0., 0.22, 0.78},
     {0., 0., 0.1}}
    >>> ma = MomentsFromMG[a, A, 5];
    >>> Print[ma];
    {4.93827160493827, 34.80707678238542, 339.49243437478106, 4335.792870444855, 68954.07073262692}
    >>> mb = MomentsFromMG[b, B, 5];
    >>> Print[mb];
    {4.93827160493827, 34.807076782385415, 339.492434374781, 4335.792870444855, 68954.07073262692}
    >>> flag = CheckDPHRepresentation[b, B];
    >>> Print[flag];
    True
    >>> Cm = SimilarityMatrix[A, B];
    >>> err1 = Norm[A.Cm-Cm.B];
    >>> err2 = Norm[a.Cm-b];
    >>> Print[Max[err1, err2]];
    3.723687702222371*^-16

    For Python/Numpy:

    >>> a = ml.matrix([[0,0,1.0]])
    >>> A = ml.matrix([[0.22, 0, 0],[0.3, 0.1, 0.55],[0.26, 0, 0.73]])
    >>> b, B = AcyclicDPHFromMG(a, A)
    >>> print(b)
    [[  1.12101e-16   9.62963e-01   3.70370e-02]]
    >>> print(B)
    [[ 0.10+0.j  0.90+0.j  0.00+0.j]
     [ 0.00+0.j  0.22+0.j  0.78+0.j]
     [ 0.00+0.j  0.00+0.j  0.73+0.j]]
    >>> ma = MomentsFromMG(a, A, 5)
    >>> print(ma)
    [4.9382716049382713, 34.807076782385423, 339.49243437478106, 4335.792870444855, 68954.070732626918]
    >>> mb = MomentsFromMG(b, B, 5)
    >>> print(mb)
    [(4.938271604938274+0j), (34.807076782385444+0j), (339.49243437478117+0j), (4335.7928704448568+0j), (68954.070732626948+0j)]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    8.23408130779e-16

