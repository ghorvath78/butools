Tools for Phase-Type Distributions (:mod:`butools.ph`)
======================================================

.. module:: butools.ph

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/ph')` 
      - in Matlab,
    * - :code:`<<BuTools`PH`` 
      - in Mathematica,
    * - :code:`from butools.ph import *` 
      - in Python/Numpy.

Phase-type and matrix-exponential distributions
-----------------------------------------------

Continuous time phase-type (PH) distributions are characterized by two parameters, the initial probability vector :math:`\alpha`
and the transient generator matrix of a continuous time (transient) Markov chain :math:`A`. The PH distribution represents
the absorption time of this transient Markov chain starting from :math:`\alpha`.
The cummulative distribution function (cdf) of PH distributions is 

.. math::
    F_{PH}(t) = 1-\alpha e^{A\,t} \mathbf{1},
    
where :math:`\mathbf{1}` is the column vector of ones.

Matrix exponential (ME) distributions are the generalizations of PH distributions. Formally, the cdf and all the related formulas for 
the statistical quantities are very similar to the ones of PH distributions

.. math::
    F_{ME}(t) = 1-\mathbf{b} e^{B\,t} \mathbf{e},
    
however, :math:`\mathbf{b},B` and :math:`\mathbf{e}` can hold general numbers, the entries do not have to be valid probabilities or
transition rates. ME distributions therefore lack the simple stochastic interpretation that PH distributions have. Vector 
:math:`\mathbf{b}` is called "starting operator", matrix :math:`B` is the "process rate operator" and vector :math:`\mathbf{e}` 
is the "summing operator" [1]_.

BuTools provides several tools for both distribution classes in the ph package. However, BuTools uses a special form of ME 
distributions, where the summing operator is a vector of ones, thus :math:`\mathbf{e}=\mathbf{1}`. This is not a restriction,
as all ME distributions defined with general summing operator can be easily transformed to this representation. Assume we
have an ME distribution in the general form with parameters :math:`\mathbf{b}, B, \mathbf{e}`. The necessary similarity transform
is obtained by calling the :code:`T = TransformToOnes(e)` procedure of the BuTools reptrans package. The parameters of the ME
distribution used by all related BuTools tools can be calculated by :math:`\mathbf{b'}=\mathbf{b}T^{-1}` and 
:math:`B'=T\,B\,T^{-1}`.

Simple statistical properties and tools
---------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MomentsFromME <butools.ph.MomentsFromME>`
      - Returns the moments of a matrix-exponential distribution.
    * - :py:func:`MomentsFromPH <butools.ph.MomentsFromPH>`
      - Returns the moments of a continuous phase-type distribution.
    * - :py:func:`CdfFromME <butools.ph.CdfFromME>`
      - Returns the cummulative distribution function of a matrix-exponential distribution.
    * - :py:func:`CdfFromPH <butools.ph.CdfFromPH>`
      - Returns the cummulative distribution function of a continuous phase-type distribution.
    * - :py:func:`PdfFromME <butools.ph.PdfFromME>`
      - Returns the probability density function of a matrix-exponential distribution.
    * - :py:func:`PdfFromPH <butools.ph.PdfFromPH>`
      - Returns the probability density function of a continuous phase-type distribution.
    * - :py:func:`IntervalPdfFromPH <butools.ph.IntervalPdfFromPH>`
      - Returns the approximate probability density function of a continuous phase-type distribution, based on the probability of falling into intervals.
    * - :py:func:`RandomPH <butools.ph.RandomPH>`
      - Returns a random phase-type distribution with a given order.
    * - :py:func:`CheckPHRepresentation <butools.ph.CheckPHRepresentation>`
      - Checks if the given vector and matrix define a valid phase-type representation.
    * - :py:func:`CheckMERepresentation <butools.ph.CheckMERepresentation>`
      - Checks if the given vector and matrix define a valid matrix-exponential representation.
    * - :py:func:`CheckMEPositiveDensity <butools.ph.CheckMEPositiveDensity>`
      - Checks if the given matrix-exponential distribution has positive density.
    * - :py:func:`SamplesFromPH <butools.ph.SamplesFromPH>`
      - Generates random samples from a phase-type distribution.
    * - :py:func:`ImageFromPH <butools.ph.ImageFromPH>`
      - Depicts the given phase-type distribution, and either displays it or saves it to file.

Inverse characterization tools
------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`APHFrom2Moments <butools.ph.APHFrom2Moments>`
      - Returns an acyclic PH which has the same 2 moments as given.
    * - :py:func:`APHFrom3Moments <butools.ph.APHFrom3Moments>`
      - Returns an acyclic PH which has the same 3 moments as given.
    * - :py:func:`PH2From3Moments <butools.ph.PH2From3Moments>`
      - Returns a PH(2) which has the same 3 moments as given.
    * - :py:func:`PH3From5Moments <butools.ph.PH3From5Moments>`
      - Returns a PH(3) which has the same 5 moments as given.
    * - :py:func:`MEFromMoments <butools.ph.MEFromMoments>`
      - Creates a matrix-exponential distribution that has the same moments as given.     
    * - :py:func:`APH2ndMomentLowerBound <butools.ph.APH2ndMomentLowerBound>`
      - Returns the lower bound of the second moment of acyclic phase-type (APH) distributions of order n.
    * - :py:func:`APH3rdMomentLowerBound <butools.ph.APH3rdMomentLowerBound>`
      - Returns the lower bound of the third moment of acyclic phase-type (APH) distributions of order n.
    * - :py:func:`APH3rdMomentUpperBound <butools.ph.APH3rdMomentUpperBound>`
      - Returns the upper bound of the third moment of acyclic phase-type (APH) distributions of order n.

Representation transformation methods
-------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CanonicalFromPH2 <butools.ph.CanonicalFromPH2>`
      - Returns the canonical form of an order-2 phase-type distribution.
    * - :py:func:`CanonicalFromPH3 <butools.ph.CanonicalFromPH3>`
      - Returns the canonical form of an order-3 phase-type distribution.
    * - :py:func:`AcyclicPHFromME <butools.ph.AcyclicPHFromME>`
      - Transforms an arbitrary matrix-exponential representation to an acyclic phase-type representation. 
    * - :py:func:`MonocyclicPHFromME <butools.ph.MonocyclicPHFromME>`
      - Transforms an arbitrary matrix-exponential representation to a Markovian monocyclic representation.
    * - :py:func:`PHFromME <butools.ph.PHFromME>`
      - Obtains a Markovian representation of a matrix exponential distribution of the same size, if possible.

Minimal representations
-----------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MEOrder <butools.ph.MEOrder>`
      - Returns the order of the ME distribution (which is not necessarily equal to the size of the representation).
    * - :py:func:`MEOrderFromMoments <butools.ph.MEOrderFromMoments>`
      - Returns the order of ME distribution that can realize the given moments.
    * - :py:func:`MinimalRepFromME <butools.ph.MinimalRepFromME>`
      - Returns the minimal representation of the given ME distribution.

.. toctree::
    :hidden:

    MomentsFromME
    MomentsFromPH
    CdfFromME
    CdfFromPH
    PdfFromME
    PdfFromPH
    IntervalPdfFromPH
    RandomPH
    CheckPHRepresentation
    CheckMERepresentation
    CheckMEPositiveDensity
    SamplesFromPH
    ImageFromPH
    APHFrom2Moments
    APHFrom3Moments
    PH2From3Moments
    PH3From5Moments
    MEFromMoments
    APH2ndMomentLowerBound
    APH3rdMomentLowerBound
    APH3rdMomentUpperBound
    AcyclicPHFromME
    CanonicalFromPH2
    CanonicalFromPH3
    MonocyclicPHFromME
    PHFromME
    MEOrder
    MEOrderFromMoments
    MinimalRepFromME
    
    
References
----------

.. [1] L. Lipsky, Queueing Theory: A Linear Algebraic Approach, MacMillan, New York, 1992.
