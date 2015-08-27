Markov Chains (:mod:`butools.mc`)
=================================

.. module:: butools.mc

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/mc')` 
      - in Matlab,
    * - :code:`<<BuTools`MC` 
      - in Mathematica,
    * - :code:`from butools.mc import *` 
      - in Python/Numpy.

Markov chains and rational processes
------------------------------------

The difference between Markov chains and rational processes is that
the entries of the generator of a rational process can be negative 
as well. There are still restrictions: the rowsum must be 0 or 1 
(in the continuous time and discrete time cases, respectively), the eigenvalues
must fall into the left half-plane in the continuous case, and in the unit cycle
in the discrete case. The dominant eigenvalue must be real.

The invariant (stationary) vector of Markov chains and rational processes
are obtained the same way. The only reason to treat them separately is 
checking the input parameters. If the global variable :code:`butools.checkInput` 
is set to :code:`true`, CTMCSolve and DTMCSolve enforce a proper Markovian generator,
while CRPSolve and DRPSolve just check the eigenvalues and the rowsums.


Properties of Markov chains and rational processes
--------------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CTMCSolve <butools.mc.CTMCSolve>`
      - Returns the steady state solution of a continuous time Markov chain
    * - :py:func:`DTMCSolve <butools.mc.DTMCSolve>`
      - Returns the steady state solution of a discrete time Markov chain
    * - :py:func:`CRPSolve <butools.mc.CRPSolve>`
      - Returns the steady state solution of a continuous time rational process
    * - :py:func:`DRPCSolve <butools.mc.DRPSolve>`
      - Returns the steady state solution of a discrete time rational process
    * - :py:func:`CheckGenerator <butools.mc.CheckGenerator>`
      - Checks if a matrix is a valid generator of a CTMC
    * - :py:func:`CheckProbMatrix <butools.mc.CheckProbMatrix>`
      - Checks if a matrix is a valid transition probability matrix of a DTMC
    * - :py:func:`CheckProbVector <butools.mc.CheckProbVector>`
      - Checks if a vector is a valid probability vector


.. toctree::
    :hidden:

    CTMCSolve
    DTMCSolve
    CRPSolve
    DRPSolve
    CheckGenerator
    CheckProbMatrix
    CheckProbVector

