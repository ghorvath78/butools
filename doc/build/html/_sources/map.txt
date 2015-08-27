Tools for Markovian Arrival Processes (:mod:`butools.map`)
==========================================================

.. module:: butools.map

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/map')` 
      - in Matlab,
    * - :code:`<<BuTools`MAP`` 
      - in Mathematica,
    * - :code:`from butools.map import *` 
      - in Python/Numpy.

Markovian arrival processes and rational arrival processes
----------------------------------------------------------

Continuous time Markovian arrival processes (MAPs) are characterized by two matrices, :math:`D_0` and :math:`D_1`. MAPs have
a continuous time Markov chain in the background with generator :math:`D=D_0+D_1`. Every time this Markov chain changes state
through a transition belonging to :math:`D_1`, an arrival event is generated. As the state of the background process is not
re-initiated after the arrival events, MAPs are capable of generating correlated arrivals.

Rational arrival processes (RAPs, also known as matrix-exponential processes, MEPs) are the generalizations of MAPs. Formally, 
all formulas for the statistical quantities are very similar to the ones of MAPs. However, :math:`D_0` and :math:`D_1` can hold 
general numbers, the entries do not have to be valid transition rates. RAPs therefore lack the simple stochastic interpretation 
that MAPs have.

Both MAPs and RAPs can be generalized to multi-type arrival processes. If there are K different arrival types, marked 
MAPs (MMAPs) and marked RAPs (MRAPs) defined by matrices :math:`D_0,\dots,D_K`  are able to describe the multi-type arrival process.

BuTools provides several tools for MAPs, RAPs and their marked variants in the maps package.

Simple statistical properties and tools
---------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MarginalDistributionFromMAP <butools.map.MarginalDistributionFromMAP>`
      - Returns the phase type marginal distribution of a Markovian arrival process.
    * - :py:func:`MarginalDistributionFromRAP <butools.map.MarginalDistributionFromRAP>`
      - Returns the matrix-exponential marginal distribution of a rational arrival process.
    * - :py:func:`MarginalDistributionFromMMAP <butools.map.MarginalDistributionFromMMAP>`
      - Returns the phase type marginal distribution of a marked Markovian arrival process.     
    * - :py:func:`MarginalDistributionFromMRAP <butools.map.MarginalDistributionFromMRAP>`
      - Returns the matrix-exponential marginal distribution of a marked rational arrival process.
    * - :py:func:`MarginalMomentsFromMAP <butools.map.MarginalMomentsFromMAP>`
      - Returns the moments of the marginal distribution of a Markovian arrival process.
    * - :py:func:`MarginalMomentsFromRAP <butools.map.MarginalMomentsFromRAP>`
      - Returns the moments of the marginal distribution of a rational arrival process.
    * - :py:func:`MarginalMomentsFromMMAP <butools.map.MarginalMomentsFromMMAP>`
      - Returns the moments of the marginal distribution of a marked Markovian arrival process.     
    * - :py:func:`MarginalMomentsFromMRAP <butools.map.MarginalMomentsFromMRAP>`
      - Returns the moments of the marginal distribution of a marked rational arrival process.
    * - :py:func:`LagCorrelationsFromMAP <butools.map.LagCorrelationsFromMAP>`
      - Returns the lag autocorrelations of a Markovian arrival process.
    * - :py:func:`LagCorrelationsFromRAP <butools.map.LagCorrelationsFromRAP>`
      - Returns the lag autocorrelations of a rational arrival process.
    * - :py:func:`LagkJointMomentsFromMAP <butools.map.LagkJointMomentsFromMAP>`
      - Returns the lag-k joint moments of a Markovian arrival process.
    * - :py:func:`LagkJointMomentsFromRAP <butools.map.LagkJointMomentsFromRAP>`
      - Returns the lag-k joint moments of a rational arrival process.
    * - :py:func:`LagkJointMomentsFromMMAP <butools.map.LagkJointMomentsFromMMAP>`
      - Returns the lag-k joint moments of a marked Markovian arrival process.
    * - :py:func:`LagkJointMomentsFromMRAP <butools.map.LagkJointMomentsFromMRAP>`
      - Returns the lag-k joint moments of a marked rational arrival process.
    * - :py:func:`RandomMAP <butools.map.RandomMAP>`
      - Returns a random Markovian arrival process with given mean value.
    * - :py:func:`RandomMMAP <butools.map.RandomMMAP>`
      - Returns a random marked Markovian arrival process with given mean value.
    * - :py:func:`CheckMAPRepresentation <butools.map.CheckMAPRepresentation>`
      - Checks if the input matrixes define a continuous time MAP.
    * - :py:func:`CheckRAPRepresentation <butools.map.CheckRAPRepresentation>`
      - Checks if the input matrixes define a continuous time RAP.
    * - :py:func:`CheckMMAPRepresentation <butools.map.CheckMMAPRepresentation>`
      - Checks if the input matrixes define a continuous time MMAP.
    * - :py:func:`CheckMRAPRepresentation <butools.map.CheckMRAPRepresentation>`
      - Checks if the input matrixes define a continuous time MRAP.
    * - :py:func:`SamplesFromMAP <butools.map.SamplesFromMAP>`
      - Generates random samples from a Markovian arrival process.
    * - :py:func:`SamplesFromMMAP <butools.map.SamplesFromMMAP>`
      - Generates random samples from a marked Markovian arrival process.
    * - :py:func:`ImageFromMAP <butools.map.ImageFromMAP>`
      - Depicts the given Markovian arrival process, and either displays it or saves it to file.
    * - :py:func:`ImageFromMMAP <butools.map.ImageFromMMAP>`
      - Depicts the given marked Markovian arrival process, and either displays it or saves it to file.
      

Inverse characterization tools
------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`RAPFromMoments <butools.map.RAPFromMoments>`
      - Creates a rational arrival process that has the same marginal and lag-1 joint moments as given.
    * - :py:func:`MRAPFromMoments <butools.map.MRAPFromMoments>`
      - Creates a marked rational arrival process that has the same marginal and lag-1 joint moments as given.
    * - :py:func:`RAPFromMomentsAndCorrelations <butools.map.RAPFromMomentsAndCorrelations>`
      - Returns a rational arrival process that has the same moments and lag autocorrelation coefficients as given.
    * - :py:func:`MAP2FromMoments <butools.map.MAP2FromMoments>`
      - Returns a MAP(2) which has the same 3 marginal moments and lag-1 autocorrelation as given.
    * - :py:func:`MAP2CorrelationBounds <butools.map.MAP2CorrelationBounds>`
      - Returns the upper and lower correlation bounds for a MAP(2) given the three marginal moments.
    * - :py:func:`MAPFromFewMomentsAndCorrelations <butools.map.MAPFromFewMomentsAndCorrelations>`
      - Returns a MAP that matches the given 2 or 3 moments and the lag-1 autocorrelation.


Representation transformation methods
-------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CanonicalFromMAP2 <butools.map.CanonicalFromMAP2>`
      - Returns the canonical form of an order-2 Markovian arrival process.
    * - :py:func:`MAPFromRAP <butools.map.MAPFromRAP>`
      - Obtains a Markovian representation of a rational arrival process of the same size, if possible.
    * - :py:func:`MMAPFromMRAP <butools.map.MMAPFromMRAP>`
      - Obtains a Markovian representation of a marked rational arrival process of the same size, if possible.
    * - :py:func:`MinimalRepFromRAP <butools.map.MinimalRepFromRAP>`
      - Returns the minimal representation of a rational arrival process.
    * - :py:func:`MinimalRepFromMRAP <butools.map.MinimalRepFromMRAP>`
      - Returns the minimal representation of a marked rational arrival process.


.. toctree::
    :hidden:

    MarginalDistributionFromMAP
    MarginalDistributionFromRAP
    MarginalDistributionFromMMAP
    MarginalDistributionFromMRAP
    MarginalMomentsFromMAP
    MarginalMomentsFromRAP
    MarginalMomentsFromMMAP
    MarginalMomentsFromMRAP
    LagCorrelationsFromMAP
    LagCorrelationsFromRAP
    LagkJointMomentsFromMAP
    LagkJointMomentsFromRAP
    LagkJointMomentsFromMMAP
    LagkJointMomentsFromMRAP
    RandomMAP
    RandomMMAP
    CheckMAPRepresentation
    CheckRAPRepresentation
    CheckMMAPRepresentation
    CheckMRAPRepresentation
    SamplesFromMAP
    SamplesFromMMAP
    ImageFromMAP
    ImageFromMMAP
    RAPFromMoments
    MRAPFromMoments
    RAPFromMomentsAndCorrelations
    MAP2FromMoments
    MAP2CorrelationBounds
    MAPFromFewMomentsAndCorrelations
    CanonicalFromMAP2
    MAPFromRAP
    MMAPFromMRAP
    MinimalRepFromRAP
    MinimalRepFromMRAP
    
