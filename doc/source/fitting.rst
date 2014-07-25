Tools for fitting measurement traces (:mod:`butools.fitting`)
=============================================================

.. currentmodule:: butools.fitting

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/fitting')` 
      - in Matlab,
    * - :code:`<<"BuTools`Fitting"` 
      - in Mathematica,
    * - :code:`from butools.fitting import *` 
      - in Python/Numpy.

Overview
--------

In the Fitting package BuTools provides functions both for fitting and
for the evaluation of the success of the fitting.

For PH and MAP fitting, currently only EM algorithms are implemented, however,
several matching procedures can also be found in the map and ph packages.

To evaluate the fitting results there are functions to compute
the distance between continuous functions given by samples (like pdf, cdf, ccdf, etc.),
and between discrtete functions (like autocorrelation function, etc.).
The efficient calculation of the log-likelihood is also supported.



Functions for fitting
---------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`PHFromTrace <butools.fitting.PHFromTrace>`
      - Performs PH distribution fitting using the EM algorithm
    * - :py:func:`MAPFromTrace <butools.fitting.MAPFromTrace>`
      - Performs MAP fitting using the EM algorithm

Functions for evaluating the results of fitting
-----------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`SquaredDifference <butools.fitting.SquaredDifference>`
      - Returns the squared difference between two vectors
    * - :py:func:`RelativeEntropy <butools.fitting.RelativeEntropy>`
      - Returns the relative entropy (aka Kullback–Leibler divergence) of two vectors
    * - :py:func:`EmpiricalSquaredDifference <butools.fitting.EmpiricalSquaredDifference>`
      - Returns the squared difference of two continuous functions given by samples and the bounds of the corresponding intervalls
    * - :py:func:`EmpiricalRelativeEntropy <butools.fitting.EmpiricalRelativeEntropy>`
      - Returns the relative entropy (aka Kullback–Leibler divergence) of two continuous functions given by samples and the bounds of the corresponding intervalls
    * - :py:func:`LikelihoodFromTrace <butools.fitting.LikelihoodFromTrace>`
      - Evaluates the log-likelihood of a trace with the given PH distribution or MAP


.. toctree::
    :hidden:

    PHFromTrace
    MAPFromTrace
    SquaredDifference
    RelativeEntropy
    EmpiricalSquaredDifference
    EmpiricalRelativeEntropy
    LikelihoodFromTrace

