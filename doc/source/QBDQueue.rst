butools.queues.QBDQueue
=======================

.. currentmodule:: butools.queues

.. np:function:: QBDQueue

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = QBDQueue(B, L, F, L0, ...)`
        * - Mathematica:
          - :code:`Ret = QBDQueue[B, L, F, L0, ...]`
        * - Python/Numpy:
          - :code:`Ret = QBDQueue(B, L, F, L0, ...)`

    Returns various performane measures of a QBD queue.

    QBD queues have a background continuous time Markov chain
    with generator Q whose the transitions can be partitioned
    into three sets: transitions accompanied by an arrival
    of a new job (F, forward), transitions accompanied by 
    the service of the current job in the server (B, 
    backward) and internal transitions (L, local). 
    Thus we have Q=B+L+F.

    Parameters
    ----------
    B : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        the service of the current job in the server
    L : matrix, shape(N,N)
        Internal transitions of the background process 
        that do not generate neither arrival nor service
    F : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        an arrival of a new job
    L0 : matrix, shape(N,N)
        Internal transitions of the background process when
        there are no jobs in the queue

    The rest of the function parameters specify the options
    and the performance measures to be computed.

    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.

The supported performance measures in this function are:

+----------------+--------------------+--------------------------------------+
| Parameter name | Input parameters   | Output                               |
+================+====================+======================================+
| "qlMoms"       | Number of moments  | The queue length moments             |
+----------------+--------------------+--------------------------------------+
| "qlDistr"      | A vector of points | The queue length distribution at     |
|                |                    | the requested points                 |
+----------------+--------------------+--------------------------------------+
| "qlDistrMG"    | None               | The vector-matrix parameters of the  |
|                |                    | matrix-geometrically distributed     |
|                |                    | queue length distribution            |
+----------------+--------------------+--------------------------------------+
| "qlDistrDPH"   | None               | The vector-matrix parameters of the  |
|                |                    | matrix-geometrically distributed     |
|                |                    | queue length distribution, converted |
|                |                    | to a discrete PH representation      |
+----------------+--------------------+--------------------------------------+
| "stMoms"       | Number of moments  | The sojourn time moments             |
+----------------+--------------------+--------------------------------------+
| "stDistr"      | A vector of points | The sojourn time distribution at the |
|                |                    | requested points (cummulative, cdf)  |
+----------------+--------------------+--------------------------------------+
| "stDistrME"    | None               | The vector-matrix parameters of the  |
|                |                    | matrix-exponentially distributed     |
|                |                    | sojourn time distribution            |
+----------------+--------------------+--------------------------------------+
| "stDistrPH"    | None               | The vector-matrix parameters of the  |
|                |                    | matrix-exponentially distributed     |
|                |                    | sojourn time distribution, converted |
|                |                    | to a continuous PH representation    |
+----------------+--------------------+--------------------------------------+
| "prec"         | The precision      | Numerical precision to check if the  |
|                |                    | input is valid and it is also used   |
|                |                    | as a stopping condition when solving |
|                |                    | the matrix-quadratic equation        |
+----------------+--------------------+--------------------------------------+

(The queue length related quantities include the customer 
in the server, and the sojourn time related quantities 
include the service times as well)

Note that "qlDistrMG" and "stDistrMG" behave much better numerically than 
"qlDistrDPH" and "stDistrPH".


