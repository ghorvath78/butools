Tools for Discrete Markovian Arrival Processes (:mod:`butools.dmap`)
====================================================================

.. module:: butools.dmap

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/dmap')` 
      - in Matlab,
    * - :code:`<<BuTools`DMAP`` 
      - in Mathematica,
    * - :code:`from butools.dmap import *` 
      - in Python/Numpy.

Discrete Markovian arrival processes and rational arrival processes
-------------------------------------------------------------------

Discrete time Markovian arrival processes (DMAPs) are characterized by two matrices, :math:`D_0` and :math:`D_1`. DMAPs have
a discrete time Markov chain in the background with transition probability matrix :math:`D=D_0+D_1`. Every time this Markov chain changes state
through a transition belonging to :math:`D_1`, an arrival event is generated. As the state of the background process is not
re-initiated after the arrival events, DMAPs are capable of generating correlated arrival times.

Discrete rational arrival processes (DRAPs) are the generalizations of DMAPs. Formally, 
all formulas for the statistical quantities are very similar to the ones of DMAPs. However, :math:`D_0` and :math:`D_1` can hold 
general numbers, the entries do not have to be valid probabilities. DRAPs therefore lack the simple stochastic interpretation 
that DMAPs have.

Both DMAPs and DRAPs can be generalized to multi-type arrival processes. If there are K different arrival types, marked 
DMAPs (DMMAPs) and marked DRAPs (DMRAPs) defined by matrices :math:`D_0,\dots,D_K`  are able to describe the multi-type arrival process.

BuTools provides several tools for DMAPs, DRAPs and their marked variants in the dmaps package.

Simple statistical properties and tools
---------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MarginalDistributionFromDMAP <butools.dmap.MarginalDistributionFromDMAP>`
      - Returns the discrete phase type distributed marginal distribution of a discrete Markovian arrival process.
    * - :py:func:`MarginalDistributionFromDRAP <butools.dmap.MarginalDistributionFromDRAP>`
      - Returns the matrix-geometrically distributed marginal of a discrete rational arrival process.
    * - :py:func:`MarginalDistributionFromDMRAP <butools.dmap.MarginalDistributionFromDMRAP>`
      - Returns the matrix-geometrically distributed marginal of a discrete marked rational arrival process.
    * - :py:func:`MarginalDistributionFromDMMAP <butools.dmap.MarginalDistributionFromDMMAP>`
      - Returns the discrete phase type distributed marginal of a discrete marked Markovian arrival process.     
    * - :py:func:`MarginalMomentsFromDMAP <butools.dmap.MarginalMomentsFromDMAP>`
      - Returns the moments of the marginal distribution of a discrete Markovian arrival process.
    * - :py:func:`MarginalMomentsFromDRAP <butools.dmap.MarginalMomentsFromDRAP>`
      - Returns the moments of the marginal distribution of a discrete rational arrival process.
    * - :py:func:`MarginalMomentsFromDMMAP <butools.dmap.MarginalMomentsFromDMMAP>`
      - Returns the moments of the marginal distribution of a discrete marked Markovian arrival process.     
    * - :py:func:`MarginalMomentsFromDMRAP <butools.dmap.MarginalMomentsFromDMRAP>`
      - Returns the moments of the marginal distribution of a discrete marked rational arrival process.
    * - :py:func:`LagCorrelationsFromDMAP <butools.dmap.LagCorrelationsFromDMAP>`
      - Returns the lag autocorrelations of a discrete Markovian arrival process.
    * - :py:func:`LagCorrelationsFromDRAP <butools.dmap.LagCorrelationsFromDRAP>`
      - Returns the lag autocorrelations of a discrete rational arrival process.
    * - :py:func:`LagkJointMomentsFromDMAP <butools.dmap.LagkJointMomentsFromDMAP>`
      - Returns the lag-k joint moments of a discrete Markovian arrival process.
    * - :py:func:`LagkJointMomentsFromDRAP <butools.dmap.LagkJointMomentsFromDRAP>`
      - Returns the lag-k joint moments of a discrete rational arrival process.
    * - :py:func:`LagkJointMomentsFromDMMAP <butools.dmap.LagkJointMomentsFromDMMAP>`
      - Returns the lag-k joint moments of a discrete marked Markovian arrival process.
    * - :py:func:`LagkJointMomentsFromDMRAP <butools.dmap.LagkJointMomentsFromDMRAP>`
      - Returns the lag-k joint moments of a discrete marked rational arrival process.
    * - :py:func:`RandomDMAP <butools.dmap.RandomDMAP>`
      - Returns a random discrete Markovian arrival process.
    * - :py:func:`RandomDMMAP <butools.dmap.RandomDMMAP>`
      - Returns a random discrete marked Markovian arrival process.
    * - :py:func:`CheckDMAPRepresentation <butools.dmap.CheckDMAPRepresentation>`
      - Checks if the input matrixes define a discrete time MAP.
    * - :py:func:`CheckDRAPRepresentation <butools.dmap.CheckDRAPRepresentation>`
      - Checks if the input matrixes define a discrete time RAP.
    * - :py:func:`CheckDMMAPRepresentation <butools.dmap.CheckDMMAPRepresentation>`
      - Checks if the input matrixes define a discrete time MMAP.
    * - :py:func:`CheckDMRAPRepresentation <butools.dmap.CheckDMRAPRepresentation>`
      - Checks if the input matrixes define a discrete time MRAP.
    * - :py:func:`SamplesFromDMAP <butools.dmap.SamplesFromDMAP>`
      - Generates random samples from a discrete Markovian arrival process.
    * - :py:func:`SamplesFromDMMAP <butools.dmap.SamplesFromDMMAP>`
      - Generates random samples from a discrete marked Markovian arrival process.
    * - :py:func:`ImageFromDMAP <butools.dmap.ImageFromDMAP>`
      - Depicts the given discrete Markovian arrival process, and either displays it or saves it to file.
    * - :py:func:`ImageFromDMMAP <butools.dmap.ImageFromDMMAP>`
      - Depicts the given discrete marked Markovian arrival process, and either displays it or saves it to file.
      

Inverse characterization tools
------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`DRAPFromMoments <butools.dmap.DRAPFromMoments>`
      - Creates a discrete rational arrival process that has the same marginal and lag-1 joint moments as given.
    * - :py:func:`DMRAPFromMoments <butools.dmap.DMRAPFromMoments>`
      - Creates a discrete marked rational arrival process that has the same marginal and lag-1 joint moments as given.
    * - :py:func:`DMAP2FromMoments <butools.dmap.DMAP2FromMoments>`
      - Returns a DMAP(2) which has the same 3 marginal moments and lag-1 autocorrelation as given.


Representation transformation methods
-------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`CanonicalFromDMAP2 <butools.dmap.CanonicalFromDMAP2>`
      - Returns the canonical form of an order-2 discrete Markovian arrival process.
    * - :py:func:`DMAPFromDRAP <butools.dmap.DMAPFromDRAP>`
      - Obtains a Markovian representation of a discrete rational arrival process of the same size, if possible.
    * - :py:func:`DMMAPFromDMRAP <butools.dmap.DMMAPFromDMRAP>`
      - Obtains a Markovian representation of a discrete marked rational arrival process of the same size, if possible.


.. toctree::
    :hidden:

    MarginalDistributionFromDMAP
    MarginalDistributionFromDRAP
    MarginalDistributionFromDMMAP
    MarginalDistributionFromDMRAP
    MarginalMomentsFromDMAP
    MarginalMomentsFromDRAP
    MarginalMomentsFromDMMAP
    MarginalMomentsFromDMRAP
    LagCorrelationsFromDMAP
    LagCorrelationsFromDRAP
    LagkJointMomentsFromDMAP
    LagkJointMomentsFromDRAP
    LagkJointMomentsFromDMMAP
    LagkJointMomentsFromDMRAP
    RandomDMAP
    RandomDMMAP
    CheckDMAPRepresentation
    CheckDRAPRepresentation
    CheckDMMAPRepresentation
    CheckDMRAPRepresentation
    SamplesFromDMAP
    SamplesFromDMMAP
    ImageFromDMAP
    ImageFromDMMAP
    DRAPFromMoments
    DMRAPFromMoments
    DMAP2FromMoments
    CanonicalFromDMAP2
    DMAPFromDRAP
    DMMAPFromDMRAP
    
