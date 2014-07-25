Tools for the Analysis of Queues (:mod:`butools.queues`)
========================================================

.. currentmodule:: butools.queues

To load this package, either start the :func:`BuToolsInit()` 
script, or execute

.. list-table::
    :widths: 50 50

    * - :code:`addpath('butools/queues')` 
      - in Matlab,
    * - :code:`<<"BuTools`Queues"` 
      - in Mathematica,
    * - :code:`from butools.queues import *` 
      - in Python/Numpy.

Type of queues supported by BuTools
-----------------------------------
The BuTools queues package supports the efficient analysis of 
queueing models that are driven by some Markovian structures.

The supported queues can be classified to two groups:

* Discrete queues, where the jobs waiting in the queue are discrete.
  In most of the supported queueing models the inter-arrival and the service times are either given by 
  Markovian arrival processes or PH distributions. 
  BuTools has a relatively large set of tools to obtain, transform and analyze MAPs and PH distributions,
  furthermore, several performance measures of these queues are also PH distributed. This kind of 
  elegant closeness property is the main motivation to restrict the focus of BuTools on these queues.
* In continuous queues the jobs are intinifesimally small, they are considered to be fluid drops.
  The rate at which fluid drops arrive and are being served are modulated by Markov chains. Again, 
  several performance measures are PH distributed, thus these queues have their well deserved places in 
  BuTools as well.

Performance measures
--------------------
Currently, the functions in the queues package calculate two performance measures:

* Queue length distribution. By our definition, the jobs in the server and in the waiting room
  are all included in the queue length.
* Sojourn time distribution. The sojourn time of the jobs, including the waiting time in the waiting room
  and the service time. This quantity is also called "system time" or "response time".

Functions for discrete queues
-----------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`QBDQueueQLD <butools.queues.QBDQueueQLD>`
      - Returns the matrix-geometric distribution of the queue length of a QBD queue.
    * - :py:func:`QBDQueueSTD <butools.queues.QBDQueueSTD>`
      - Returns the matrix-exponential distribution of the sojourn time of the jobs in a QBD queue.
    * - :py:func:`MAPMAP1QLD <butools.queues.MAPMAP1QLD>`
      - Returns the matrix-geometric distribution of the queue length of a MAP/MAP/1 queue.
    * - :py:func:`MAPMAP1STD <butools.queues.MAPMAP1STD>`
      - Returns the matrix-exponential distribution of the sojourn time of the jobs in a MAP/MAP/1 queue.
    * - :py:func:`MAPPH1QLD <butools.queues.MAPPH1QLD>`
      - Returns the matrix-geometric distribution of the queue length of a MAP/PH/1 queue.
    * - :py:func:`MAPPH1STD <butools.queues.MAPPH1STD>`
      - Returns the matrix-exponential distribution of the sojourn time of the jobs in a MAP/PH/1 queue.
    * - :py:func:`PHPH1QLD <butools.queues.PHPH1QLD>`
      - Returns the matrix-geometric distribution of the queue length of a PH/PH/1 queue.
    * - :py:func:`PHPH1STD <butools.queues.PHPH1STD>`
      - Returns the matrix-exponential distribution of the sojourn time of the jobs in a PH/PH/1 queue.

Functions for continuous queues
-------------------------------

.. list-table::
    :widths: 25 150

    * - :py:func:`FluidQueueQLD <butools.queues.FluidQueueQLD>`
      - Returns the matrix-exponential distribution of the fluid level in a fluid queue.
    * - :py:func:`FluidQueueSTD <butools.queues.FluidQueueSTD>`
      - Returns the matrix-exponential distribution of the sojourn time of fluid drops in a fluid queue.
    * - :py:func:`FluFluQLD <butools.queues.FluFluQLD>`
      - Returns the matrix-exponential distribution of the fluid level in a fluid queue, where the input and output processes are independent.
    * - :py:func:`FluFluSTD <butools.queues.FluFluSTD>`
      - Returns the matrix-exponential distribution of the sojourn time of fluid drops in a fluid queue, where the input and output processes are independent.

.. toctree::
    :hidden:

    QBDQueueQLD
    QBDQueueSTD
    MAPMAP1QLD
    MAPMAP1STD
    MAPPH1QLD
    MAPPH1STD
    PHPH1QLD
    PHPH1STD
    FluidQueueQLD
    FluidQueueSTD
    FluFluQLD
    FluFluSTD
 

