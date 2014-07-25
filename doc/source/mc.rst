Markov Chains (:mod:`butools.mc`)
=================================

.. currentmodule:: butools.mc

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/mc')` 
      - in Matlab,
    * - :code:`<<"BuTools`MC"` 
      - in Mathematica,
    * - :code:`from butools.mc import *` 
      - in Python/Numpy.


Properties of Markov chains and rational processes
--------------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CTMCSolve <butools.mc.CTMCSolve>`
      - Steady state solution of a continuous time Markov chain
    * - :py:func:`DTMCSolve <butools.mc.DTMCSolve>`
      - Steady state solution of a discrete time Markov chain
    * - :py:func:`CRPSolve <butools.mc.CRPSolve>`
      - Steady state solution of a continuous time rational process
    * - :py:func:`DRPCSolve <butools.mc.DRPSolve>`
      - Steady state solution of a discrete time rational process
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

