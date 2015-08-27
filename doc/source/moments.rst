Moment expressions (:mod:`butools.moments`)
===========================================

.. module:: butools.moments

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/moments')` 
      - in Matlab,
    * - :code:`<<BuTools`Moments`` 
      - in Mathematica,
    * - :code:`from butools.moments import *` 
      - in Python/Numpy.


Moment expressions
------------------

Several researchers working of fitting/matching procedures for phase-type distributions
and Markovian arrival processes have introduced moment expressions to make the 
arising expressions simpler. Instead of the classical (raw) moments, factorial, reduced, 
normalized, and Hankel moments have been defined in various research papers. BuTools has 
tools available to convert among them.

Conversion routines between various moment expressions
------------------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`NormMomsFromMoms <butools.moments.NormMomsFromMoms>`
      - Returns the normalized moments given the raw moments
    * - :py:func:`MomsFromNormMoms <butools.moments.MomsFromNormMoms>`
      - Returns the raw moments given the normalized moments
    * - :py:func:`ReducedMomsFromMoms <butools.moments.ReducedMomsFromMoms>`
      - Returns the reduced moments given the raw moments
    * - :py:func:`MomsFromReducedMoms <butools.moments.MomsFromReducedMoms>`
      - Returns the raw moments given the reduced moments
    * - :py:func:`HankelMomsFromMoms <butools.moments.HankelMomsFromMoms>`
      - Returns the Hankel moments given the raw moments
    * - :py:func:`MomsFromHankelMoms <butools.moments.MomsFromHankelMoms>`
      - Returns the raw moments given the Hankel moments
    * - :py:func:`FactorialMomsFromMoms <butools.moments.FactorialMomsFromMoms>`
      - Returns the factorial moments given the raw moments
    * - :py:func:`MomsFromFactorialMoms <butools.moments.MomsFromFactorialMoms>`
      - Returns the raw moments given the factorial moments
    * - :py:func:`JFactorialMomsFromJMoms <butools.moments.JFactorialMomsFromJMoms>`
      - Returns the joint factorial moments given the joint raw moments
    * - :py:func:`JMomsFromJFactorialMoms <butools.moments.JMomsFromJFactorialMoms>`
      - Returns the joint raw moments given the joint factorial moments
    * - :py:func:`CheckMoments <butools.moments.CheckMoments>`
      - Checks if the given moment sequence belongs to a distribution with support (0,inf)


.. toctree::
    :hidden:

    NormMomsFromMoms
    MomsFromNormMoms
    ReducedMomsFromMoms
    MomsFromReducedMoms
    HankelMomsFromMoms
    MomsFromHankelMoms
    FactorialMomsFromMoms
    MomsFromFactorialMoms
    JFactorialMomsFromJMoms
    JMomsFromJFactorialMoms
    CheckMoments

