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

BuTools is portable, no installation is needed.
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

Global variables
----------------

There are three global variables used by BuTools:

+-------------------------------+--------------------------------+--------------------------------+---------------+
| Name in MATLAB                | Name in Mathematica            | Name in Python                 | Default value |
+===============================+================================+================================+===============+
| :code:`BuToolsVerbose`        | :code:`BuTools`Verbose`        | :data:`butools.verbose`        | :code:`False` |
+-------------------------------+--------------------------------+--------------------------------+---------------+
| :code:`BuToolsCheckInput`     | :code:`BuTools`CheckInput`     | :data:`butools.checkInput`     | :code:`True`  |
+-------------------------------+--------------------------------+--------------------------------+---------------+
| :code:`BuToolsCheckPrecision` | :code:`BuTools`CheckPrecision` | :data:`butools.checkPrecision` | :code:`1e-12` |
+-------------------------------+--------------------------------+--------------------------------+---------------+

.. py:data:: verbose

    Setting :data:`verbose` to :code:`True` allows the functions to print as many useful messages to the output
    console as possible. Turning it off avoids bloating the console. The default value is :code:`False`, but for 
    the examples of the reference documentation we have set it to :code:`True`.

.. py:data:: checkInput

    If :data:`checkInput` is set to :code:`True`, the functions of BuTools perform as many checks on the input 
    parameters as possible. This can be very useful to recognize typos as soon as possible, but can be a waste
    of computational effort in case of a computationally demanding application.

.. py:data:: checkPrecision

    This numeric value serves as the tolerance when the validity of the input parameters are checked.


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

