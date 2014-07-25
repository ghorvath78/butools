Tools for the Analysis of Measurement Traces (:mod:`butools.trace`)
===================================================================

.. currentmodule:: butools.trace

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/trace')` 
      - in Matlab,
    * - :code:`<<"BuTools`Trace"` 
      - in Mathematica,
    * - :code:`from butools.trace import *` 
      - in Python/Numpy.


Tools for the analysis of measurement traces
--------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CdfFromTrace <butools.trace.CdfFromTrace>`
      - Returns the empirical distribution function of the trace
    * - :py:func:`CdfFromWeightedTrace <butools.trace.CdfFromWeightedTrace>`
      - Returns the empirical distribution function of a trace consisting of weighted data
    * - :py:func:`PdfFromTrace <butools.trace.PdfFromTrace>`
      - Returns the empirical density function of a trace
    * - :py:func:`PdfFromWeightedTrace <butools.trace.PdfFromWeightedTrace>`
      - Returns the empirical density function of a trace consisting of weighted data
    * - :py:func:`MarginalMomentsFromTrace <butools.trace.MarginalMomentsFromTrace>`
      - Returns the marginal moments of a trace
    * - :py:func:`MarginalMomentsFromWeightedTrace <butools.trace.MarginalMomentsFromWeightedTrace>`
      - Returns the marginal moments of a trace consisting of weighted data
    * - :py:func:`LagCorrelationsFromTrace <butools.trace.LagCorrelationsFromTrace>`
      - Returns the lag-k autocorrelation of a trace
    * - :py:func:`LagkJointMomentsFromTrace <butools.trace.LagkJointMomentsFromTrace>`
      - Returns the lag-k joint moments of a trace
    * - :py:func:`IATimesFromCummulative <butools.trace.IATimesFromCummulative>`
      - Returns inter-arrival times from cummulative a trace.


.. toctree::
    :hidden:

    CdfFromTrace
    CdfFromWeightedTrace
    PdfFromTrace
    PdfFromWeightedTrace
    MarginalMomentsFromTrace
    MarginalMomentsFromWeightedTrace
    LagCorrelationsFromTrace
    LagkJointMomentsFromTrace
    IATimesFromCummulative

