.. _prep_and_run:

*************************************************
Setting Up and Running the Long Baseline Pipeline
*************************************************

.. _preparation:

==================
Preparing the data
==================

.. note::

    Processing of interleaved datasets is not currently supported.

The **long baseline pipeline** requires **LOFAR HBA** pre-processed data and the solution set ``solutions.h5`` from a corresponding `prefactor`_ v3.0 run. Instructions how to run **prefactor** on your data, please look at the `prefactor documentation`_.
Data can be typically obtained from the LOFAR Long-Term Archive at https://lta.lofar.eu. All
input measurement-sets for one pipeline run need to be in the same directory.
  
.. _parset:

===============================================
Preparing the configuration and the parset file
===============================================
The **long baseline pipeline** uses the ``genericpipeline`` framwork. See the **prefafactor** `documentation`_ on how to modify the ``pipeline.cfg`` and the corresponding parset files before you start the pipeline to suit your needs.


.. _running:

====================
Running the pipeline
====================

Once all parameters are set, the pipeline can be run as, for example,::

   genericpipeline.py -c pipeline.cfg LB-Delay-Calibration.parset
   
.. _help:

============
Getting help
============

The **long baseline pipeline** is a continuously maintained and developed software package.
If you need help or want to report any issues or bugs, follow these links:

- `Long Baseline Pipeline GitHub issues`_

.. - Frequently Asked Questions (`FAQ`_)


.. _Long Baseline Pipeline GitHub issues: https://github.com/lmorabit/long_baseline_pipeline/issues
.. .. _FAQ: https://github.com/lofar-astron/prefactor/wiki/Documentation%3A-Faq  
.. _prefactor: https://github.com/lofar-astron/prefactor
.. _prefactor documentation: https://www.astron.nl/citt/prefactor/
.. _documentation: file:///media/quasarfix/media/cep3/prefactor/docs/build/html/parset.html
