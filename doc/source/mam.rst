Matrix-Analytic Methods (:mod:`butools.mam`)
============================================

.. module:: butools.mam

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/mam')` 
      - in Matlab,
    * - :code:`<<"BuTools`MAM"` 
      - in Mathematica,
    * - :code:`from butools.mam import *` 
      - in Python/Numpy.

Under Matlab, the SMCSolver package [1]_ (both for QBD and M/G/1,G/M/1 type Markov chains) must be in path as well.

Quasi Birth-Death Processes
---------------------------

Quasi Birth-Death Processes (QBDs, [2]_) are discrete time or continuous time Markov chains.
The state space is two dimensional, it is composed by the levels and the phases.
QBDs have a block-tridiagonal generator. 

.. math::
    Q=\begin{bmatrix}L_0 & F & &\\ B & L & F & \\ & B & L & F \\ & & \ddots & \ddots & \ddots \end{bmatrix}

The matrix describing the transitions to the 
previous level is denoted by :math:`B` (backward), the ones to the next level by :math:`F`
(forward), finally, transitions that are not accompanied by the change of the level
are held by matrix :math:`L` (local). 

The QBDs have a matrix-geometrically distributed stationary distribution,

.. math::
    \pi_k=\pi_0 R^k,

where the :math:`i` th entry of vector :math:`\pi_k` is the stationary probability that the
level is :math:`k` and the phase is :math:`i`.

Matrix :math:`R` is the solution of a matrix-quadratic equation.

BuTools has functions to determine both :math:`\pi_0` and :math:`R`, thus all
ingredients to easily calculate the stationary distribution of QBDs. Several functions
are just wrappers to SMCSolver, so do not forget to give credit to [1]_ as well if you are using it.


M/G/1 and G/M/1 Type Markov Chains
----------------------------------

M/G/1 and G/M/1 type Markov chains are more general than QBDs. They both have a block-triangular
generator. For M/G/1 the structure of the generator is

.. math::
    Q=\begin{bmatrix}B_0 & B_1 & B_2 & B_3 & \dots \\ A_0 & A_1 & A_2 & A_3 & \dots \\ & A_0 & A_1 & A_2 & \dots \\ & & \ddots & \ddots & \ddots \end{bmatrix}

while for G/M/1 type it is

.. math::
    Q=\begin{bmatrix}B_0 & A_0  \\ B_1 & A_1 & A_0   \\ B_2 & A_2 & A_1 & A_0  \\ \vdots & \ddots & \ddots & \ddots  & \ddots \end{bmatrix}

BuTools is able to calculate the stationary solution of such Markov chains. The
corresponding functions are mostly wrappers of SMCTools, please cite [1]_ if you are 
using them.

Markovian Fluid Models
----------------------

Markovian fluid models [3]_ are characterized by a continuous time Markov chain with generator 
:math:`Q`, and a diagonal matrix :math:`R`. The :math:`i` th entry of the diagonal of :math:`R`
is the rate at which the fluid arrives to the fluid storage when the background Markov
chain is in state :math:`i` . Fluid rates can be positive (the fluid level increases), negative 
(the fluid level decreases) and can be zero as well (while being in these states the fluid level of the 
storage remains the same).

There is a class of Markovian fluid models, called *canonical* fluid models, which are 
much easier to hande technically. In canonical fluid models the fluid rates can be either
1 or -1, zero rates or rates other than these two are not allowed.

The steady state solution of Markovian fluid models (both canonical and general ones) is 
matrix-exponential. The probability that the fluid level is exactly 0 while being in various
states of the background process are denoted by vector :math:`p_m`, and the density vector
corresponding to fluid level is :math:`x` is expressed by

.. math::
    \pi(x)=\alpha e^{Kx} C,

BuTools provides functions to obtain the initial vector :math:`\alpha`, matrix :math:`K` and
closing matrix :math:`C` both for canonical and general fluid models.


Functions to solve QBDs
-----------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`QBDFundamentalMatrices <butools.mam.QBDFundamentalMatrices>`
      - Returns any combination of the R, G, and U matrices of a discrete or continuous time QBD
    * - :py:func:`QBDSolve <butools.mam.QBDSolve>`
      - Returns the parameters of the matrix-geometrically distributed stationary distribution of a QBD
    * - :py:func:`QBDStationaryDistr <butools.mam.QBDStationaryDistr>`
      - Returns the stationary distribution of a QBD up to a given level

Functions to solve M/G/1 and G/M/1 type Markov chains
-----------------------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`MG1FundamentalMatrix <butools.mam.MG1FundamentalMatrix>`
      - Returns matrix G, the fundamental matrix of an M/G/1 type Markov chain
    * - :py:func:`MG1StationaryDistr <butools.mam.MG1StationaryDistr>`
      - Returns the stationary distribution of an M/G/1 type Markov chain up to a given level
    * - :py:func:`GM1FundamentalMatrix <butools.mam.GM1FundamentalMatrix>`
      - Returns matrix R, the fundamental matrix of a G/M/1 type Markov chain
    * - :py:func:`GM1StationaryDistr <butools.mam.GM1StationaryDistr>`
      - Returns the stationary distribution of a G/M/1 type Markov chain up to a given level


Functions to solve Markovian fluid models
-----------------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`FluidFundamentalMatrices <butools.mam.FluidFundamentalMatrices>`
      - Returns any combination of the Psi, K, and U matrices of a canonical Markovian fluid model
    * - :py:func:`FluidSolve <butools.mam.FluidSolve>`
      - Returns the parameters of the matrix-exponentially distributed stationary distribution of a canonical Markovian fluid model
    * - :py:func:`GeneralFluidSolve <butools.mam.GeneralFluidSolve>`
      - Returns the parameters of the matrix-exponentially distributed stationary distribution of a general Markovian fluid model
    * - :py:func:`FluidStationaryDistr <butools.mam.FluidStationaryDistr>`
      - Returns the cummulative distribution of the stationary fluid level of a Markovian fluid model

Refereces
---------
.. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt, B. (2006, October). Structured Markov chains solver: software tools. 
       In Proceeding from the 2006 workshop on Tools for solving structured Markov chains (p. 14). ACM.
.. [2] Latouche, Guy, Vaidyanathan Ramaswami. "Introduction to matrix geometric methods in stochastic modeling." ASA-SIAM 
       Series on Statistics and Applied Probability. SIAM, Philadelphia PA (1999).
.. [3] A da Silva Soares. "Fluid Queues." PhD dissertation, de l'Université libre de Bruxelles, 2005.


.. toctree::
    :hidden:

    QBDFundamentalMatrices
    QBDSolve
    QBDStationaryDistr
    MG1FundamentalMatrix
    MG1StationaryDistr
    GM1FundamentalMatrix
    GM1StationaryDistr
    FluidFundamentalMatrices
    FluidSolve
    GeneralFluidSolve
    FluidStationaryDistr    
    
