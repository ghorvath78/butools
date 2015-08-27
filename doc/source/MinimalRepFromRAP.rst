butools.map.MinimalRepFromRAP
=============================

.. currentmodule:: butools.map

.. np:function:: MinimalRepFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = MinimalRepFromRAP(H0, H1, how, precision)`
        * - Mathematica:
          - :code:`{D0, D1} = MinimalRepFromRAP[H0, H1, how, precision]`
        * - Python/Numpy:
          - :code:`D0, D1 = MinimalRepFromRAP(H0, H1, how, precision)`

    Returns the minimal representation of a rational arrival 
    process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    how : {"obs", "cont", "obscont"}, optional      
        Determines how the representation is minimized. "cont" 
        means controllability, "obs" means observability, 
        "obscont" means that the rational arrival process is 
        minimized in both respects. The default value is 
        "obscont"
    precision : double, optional
       Precision used by the Staircase algorithm. The default 
       value is 1e-12.
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the minimal representation
    D1 : matrix, shape (M,M)
        The D1 matrix of the minimal representation

    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.

    Examples
    ========
    For Matlab:

    >>> D0 = [-5., 1., 0; 3., -3., 0; 1., 1., -5.];
    >>> D1 = [0, 0, 4.; 0, 0, 0; 1., 1., 1.];
    >>> [H0, H1] = MinimalRepFromRAP(D0, D1, 'cont');
    >>> disp(H0);
        -5     1     0
         3    -3     0
         1     1    -5
    >>> disp(H1);
         0     0     4
         0     0     0
         1     1     1
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
       3.4773e-15
    >>> [H0, H1] = MinimalRepFromRAP(D0, D1, 'obs');
    >>> disp(H0);
          -4.4074       1.6931
          0.84259      -2.5926
    >>> disp(H1);
            2.037      0.67725
            2.787       -1.037
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
        8.806e-15
    >>> [H0, H1] = MinimalRepFromRAP(D0, D1, 'obscont');
    >>> disp(H0);
          -4.4074       1.6931
          0.84259      -2.5926
    >>> disp(H1);
            2.037      0.67725
            2.787       -1.037
    >>> Cm = SimilarityMatrix(H0, D0);
    >>> err1 = norm(H0*Cm-Cm*D0);
    >>> err2 = norm(H1*Cm-Cm*D1);
    >>> disp(max(err1, err2));
        8.806e-15

    For Mathematica:

    >>> D0 = {{-5., 1., 0},{3., -3., 0},{1., 1., -5.}};
    >>> D1 = {{0, 0, 4.},{0, 0, 0},{1., 1., 1.}};
    >>> {H0, H1} = MinimalRepFromRAP[D0, D1, "cont"];
    >>> Print[H0];
    {{-5., 1., 0.},
     {3., -3., 0.},
     {1., 1., -5.}}
    >>> Print[H1];
    {{0., 0., 4.},
     {0., 0., 0.},
     {1., 1., 1.}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    4.286666800451802*^-15
    >>> {H0, H1} = MinimalRepFromRAP[D0, D1, "obs"];
    >>> Print[H0];
    {{-4.407407407407407, 1.6931216931216935},
     {0.8425925925925922, -2.592592592592593}}
    >>> Print[H1];
    {{2.0370370370370363, 0.6772486772486779},
     {2.7870370370370363, -1.0370370370370368}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    8.396702586553632*^-15
    >>> {H0, H1} = MinimalRepFromRAP[D0, D1, "obscont"];
    >>> Print[H0];
    {{-4.407407407407407, 1.6931216931216935},
     {0.8425925925925922, -2.592592592592593}}
    >>> Print[H1];
    {{2.0370370370370363, 0.6772486772486779},
     {2.7870370370370363, -1.0370370370370368}}
    >>> Cm = SimilarityMatrix[H0, D0];
    >>> err1 = Norm[H0.Cm-Cm.D0];
    >>> err2 = Norm[H1.Cm-Cm.D1];
    >>> Print[Max[err1, err2]];
    8.396702586553632*^-15

    For Python/Numpy:

    >>> D0 = ml.matrix([[-5., 1., 0],[3., -3., 0],[1., 1., -5.]])
    >>> D1 = ml.matrix([[0, 0, 4.],[0, 0, 0],[1., 1., 1.]])
    >>> H0, H1 = MinimalRepFromRAP(D0, D1, "cont")
    >>> print(H0)
    [[-5.  1.  0.]
     [ 3. -3.  0.]
     [ 1.  1. -5.]]
    >>> print(H1)
    [[ 0.  0.  4.]
     [ 0.  0.  0.]
     [ 1.  1.  1.]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    1.82603380554e-15
    >>> H0, H1 = MinimalRepFromRAP(D0, D1, "obs")
    >>> print(H0)
    [[-4.40741  1.69312]
     [ 0.84259 -2.59259]]
    >>> print(H1)
    [[ 2.03704  0.67725]
     [ 2.78704 -1.03704]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    7.63152157294e-16
    >>> H0, H1 = MinimalRepFromRAP(D0, D1, "obscont")
    >>> print(H0)
    [[-4.40741  1.69312]
     [ 0.84259 -2.59259]]
    >>> print(H1)
    [[ 2.03704  0.67725]
     [ 2.78704 -1.03704]]
    >>> Cm = SimilarityMatrix(H0, D0)
    >>> err1 = la.norm(H0*Cm-Cm*D0)
    >>> err2 = la.norm(H1*Cm-Cm*D1)
    >>> print(np.max(err1, err2))
    7.63152157294e-16

