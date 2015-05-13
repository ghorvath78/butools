.. _reference:

#################
BuTools Reference
#################

.. py:module:: butools

About BuTools
=============

BuTools is a large and growing set of functions related to Markovian performance
and dependability analyis.

Our main motivation was to make our life easier by collecting the most important
tools which we use almost every day. The multi-platform nature of BuTools and the 
standardized call interfaces help us in the collaborative work in our research group.
Several functions do something trivial, but there
are some gems as well; procedures, which have a hard to understand background and 
difficult implementation.

We hope that BuTools we be useful not only to us but to a broader community as well.


Using BuTools
=============

Supported frameworks
--------------------

At the moment BuTools supports the following three frameworks:

* *MATLAB*: a widely used industry-standard framework with high license fees
* *NumPy/Python*: a free (open source) alternative to MATLAB, with comparable performance. The popularity of python based scientific tools is rapidly growing, a huge and helpful community is available on the internet.
* *Mathematica*: This is an expensive alternative again, but it supports arbitrary precision arithmetic and symbolic computations, that can frequently be useful.


Loading BuTools
---------------

The packages of BuTools can be loaded individually, but there are convenience functions available to load everything as well.
If BuTools is located in directory <BTDir>, all BuTools packages can be loaded by

.. list-table::
    :widths: 50 50

    * - :code:`run('<BTDir>/Matlab/BuToolsInit.m')` 
      - in Matlab,
    * - :code:`%run "<BTDir>/Python/BuToolsInit"` 
      - in an IPython console,
    * - :code:`AppendTo[$Path,"<BTDir>/Mathematica"]; <<BuTools`` 
      - in Mathematica.

In case of Matlab and IPython no other preparations and installation steps are required, BuToolsInit adjusts the Path variable, too.

The :func:`BuToolsInit()` function has two optional arguments as well, by which it is possible to initialize the BuToolsVerbose,
BuToolsCheckInput and BuToolsCheckPrecision global variables.


Packages
========

.. toctree::
   :maxdepth: 1

   utils
   mc
   moments
   reptrans
   ph
   dph
   map
   dmap
   trace
   fitting
   mam
   queues


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

