Tools for Discrete Phase-Type Distributions (:mod:`butools.dph`)
================================================================

.. currentmodule:: butools.dph

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/dph')` 
      - in Matlab,
    * - :code:`<<"BuTools`DPH"` 
      - in Mathematica,
    * - :code:`from butools.dph import *` 
      - in Python/Numpy.

Discrete phase-type and matrix-geometric distributions
------------------------------------------------------

Discrete time phase-type (DPH) distributions are characterized by two parameters, the initial probability vector :math:`\alpha`
and the transition probability matrix of a transient discrete time Markov chain :math:`\alpha`. The DPH distribution represents
the absorption time of this transient Markov chain starting from :math:`\alpha`.
The cummulative distribution function (cdf) of DPH distributions is 

.. math::
    F_{PH}(k) = P(\mathcal{X}_{PH}\leq k) = 1-\alpha A^k \mathbf{1},
    
where :math:`\mathbf{1}` is the column vector of ones.

Matrix geometric (MG) distributions are the generalizations of DPH distributions. Formally, the cdf and all the related formulas for 
the statistical quantities are very similar to the ones of DPH distributions

.. math::
    F_{MG}(t) = P(\mathcal{X}_{MG}\leq k) = 1-\mathbf{b} B^k \mathbf{e},
    
however, :math:`\mathbf{b},B` and :math:`\mathbf{e}` can hold general numbers, the entries do not have to be valid probabilities. 
MG distributions therefore lack the simple stochastic interpretation that DPH distributions have. Vector :math:`\mathbf{b}` is 
called "starting operator", matrix :math:`B` is the "process rate operator" and vector :math:`\mathbf{e}` is the "summing operator".

BuTools provides several tools for both distribution classes in the dph package. However, BuTools uses a special form of MG 
distributions, where the summing operator is a vector of ones, thus :math:`\mathbf{e}=\mathbf{1}`. This is not a restriction,
as all MG distributions defined with general summing operator can be easily transformed to this representation. Assume we
have an MG distribution in the general form with parameters :math:`\mathbf{b}, B, \mathbf{e}`. The necessary similarity transform
is obtained by calling the :code:`T = TransformToOnes(e)` procedure of the BuTools reptrans package. The parameters of the MG
distribution used by all related BuTools tools can be calculated by :math:`\mathbf{b'}=\mathbf{b}T^{-1}` and 
:math:`B'=T\,B\,T^{-1}`.

Simple statistical properties and tools
---------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MomentsFromMG <butools.dph.MomentsFromMG>`
      - Returns the moments of a matrix geometric distribution.
    * - :py:func:`MomentsFromDPH <butools.dph.MomentsFromDPH>`
      - Returns the moments of a discrete phase-type distribution.
    * - :py:func:`CdfFromMG <butools.dph.CdfFromMG>`
      - Returns the cummulative distribution function of a matrix-geometric distribution.
    * - :py:func:`CdfFromDPH <butools.dph.CdfFromDPH>`
      - Returns the cummulative distribution function of a discrete phase-type distribution.
    * - :py:func:`PmfFromMG <butools.dph.PmfFromMG>`
      - Returns the probability mass function of a matrix-geometric distribution.
    * - :py:func:`PmfFromDPH <butools.dph.PmfFromDPH>`
      - Returns the probability mass function of a discrete phase-type distribution.
    * - :py:func:`RandomDPH <butools.dph.RandomDPH>`
      - Returns a random discrete phase-type distribution with a given mean value.
    * - :py:func:`CheckDPHRepresentation <butools.dph.CheckDPHRepresentation>`
      - Checks if the given vector and matrix define a valid discrete phase-type representation.
    * - :py:func:`CheckMGRepresentation <butools.dph.CheckMGRepresentation>`
      - Checks if the given vector and matrix define a valid matrix-geometric representation.
    * - :py:func:`SamplesFromDPH <butools.dph.SamplesFromDPH>`
      - Generates random samples from a discrete phase-type distribution.
    * - :py:func:`ImageFromDPH <butools.dph.ImageFromDPH>`
      - Depicts the given discrete phase-type distribution, and either displays it or saves it to file.
      

Inverse characterization tools
------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`DPH2From3Moments <butools.dph.DPH2From3Moments>`
      - Returns an order-2 discrete phase-type distribution which has the same 3 moments as given.
    * - :py:func:`DPH3From5Moments <butools.dph.DPH3From5Moments>`
      - Returns an order-3 discrete phase-type distribution which has the same 5 moments as given.
    * - :py:func:`MGFromMoments <butools.dph.MGFromMoments>`
      - Creates a matrix-geometric distribution that has the same moments as given.

Representation transformation methods
-------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`AcyclicDPHFromMG <butools.dph.AcyclicDPHFromMG>`
      - Transforms a matrix-geometric representation to an acyclic DPH representation of the same size, if possible.
    * - :py:func:`CanonicalFromDPH2 <butools.dph.CanonicalFromDPH2>`
      - Returns the canonical form of an order-2 discrete phase-type distribution.
    * - :py:func:`CanonicalFromDPH3 <butools.dph.CanonicalFromDPH3>`
      - Returns the canonical form of an order-3 discrete phase-type distribution.
    * - :py:func:`DPHFromMG <butools.dph.DPHFromMG>`
      - Obtains a Markovian representation of a matrix geometric distribution of the same size, if possible.

.. toctree::
    :hidden:

    MomentsFromMG
    MomentsFromDPH
    CdfFromMG
    CdfFromDPH
    PmfFromMG
    PmfFromDPH
    RandomDPH
    CheckDPHRepresentation
    CheckMGRepresentation
    SamplesFromDPH
    ImageFromDPH
    DPH2From3Moments
    DPH3From5Moments
    MGFromMoments
    AcyclicDPHFromMG
    CanonicalFromDPH2
    CanonicalFromDPH3
    DPHFromMG
    
