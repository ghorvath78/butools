butools.dmap.LagkJointMomentsFromDMRAP
======================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDMRAP(H, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDMRAP[H, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDMRAP(H, K, L, prec)`

    Returns the lag-L joint moments of a discrete marked 
    rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to check
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14

    Returns
    -------
    Nm : list/cell of matrices of shape(K+1,K+1), length(L)
        The matrices containing the lag-L joint moments,
        starting from moment 0.

    Examples
    ========
    For Matlab:

    >>> H0 = [0.15, 0.2, 0.18; -0.23, 0.17, 0.22; 0.19, 0.15, 0.16];
    >>> H1 = [0.01, 0.08, 0.16; 0.02, 0.2, 0.07; 0.02, 0.15, 0.17];
    >>> H2 = [0.14, 0.07, 0.01; 0.19, 0.02, 0.34; 0.06, 0.1, 0];
    >>> Nm = LagkJointMomentsFromDMRAP({H0,H1,H2},3,2);
    >>> disp(Nm{1});
          0.48798      0.78047       1.6785       4.9029
          0.77458       1.2395       2.6673       7.7945
           1.6539       2.6481       5.7016       16.669
           4.8092       7.7033       16.593       48.526
    >>> disp(Nm{2});
          0.51202      0.81429       1.7401       5.0566
          0.82019       1.3036       2.7837       8.0853
           1.7647       2.8029       5.9814       17.365
           5.1503        8.177       17.442       50.619

    For Mathematica:

    >>> H0 = {{0.15, 0.2, 0.18},{-0.23, 0.17, 0.22},{0.19, 0.15, 0.16}};
    >>> H1 = {{0.01, 0.08, 0.16},{0.02, 0.2, 0.07},{0.02, 0.15, 0.17}};
    >>> H2 = {{0.14, 0.07, 0.01},{0.19, 0.02, 0.34},{0.06, 0.1, 0}};
    >>> Nm = LagkJointMomentsFromDMRAP[{H0,H1,H2},3,2];
    >>> Print[Nm[[1]]];
    {{0.4879805825563548, 0.7804739995572898, 1.6784777723131161, 4.902860723736972},
     {0.7745778986519039, 1.2395315649453407, 2.6672632669482894, 7.794493672457849},
     {1.653860986909174, 2.648092315164064, 5.701611221347518, 16.66906590441337},
     {4.80922115249553, 7.70333279203901, 16.592926073643227, 48.52561424626535}}
    >>> Print[Nm[[2]]];
    {{0.5120194174436454, 0.8142942617124118, 1.7400542916801724, 5.0566248687873845},
     {0.820190362617798, 1.3035610054097582, 2.7836656770888917, 8.085263953193856},
     {1.764671077084115, 2.8029046012949106, 5.98142492574121, 17.364636197177674},
     {5.1502644400288276, 8.17696708018504, 17.441960447818936, 50.61882961114629}}

    For Python/Numpy:

    >>> H0 = ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1 = ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2 = ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    >>> Nm = LagkJointMomentsFromDMRAP([H0,H1,H2],3,2)
    >>> print(Nm[0])
    [[  0.48798   0.78047   1.67848   4.90286]
     [  0.77458   1.23953   2.66726   7.79449]
     [  1.65386   2.64809   5.70161  16.66907]
     [  4.80922   7.70333  16.59293  48.52561]]
    >>> print(Nm[1])
    [[  0.51202   0.81429   1.74005   5.05662]
     [  0.82019   1.30356   2.78367   8.08526]
     [  1.76467   2.8029    5.98142  17.36464]
     [  5.15026   8.17697  17.44196  50.61883]]

