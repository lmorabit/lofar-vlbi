.. long_baseline_pipeline documentation master file, created by
   sphinx-quickstart on Mon Nov 11 17:22:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

long_baseline_pipeline: a calibration pipeline for LOFAR's international stations
=================================================================================
The LOFAR **long baseline pipeline** analyses LOFAR data incorporating all international stations. It looks and solves for delays for potential long-baseline calibrators within the field-of-view of the pointing. The pipeline provides solution tables (in h5parm format) and data sets of every potential in-field calibrator. You can use these solutions in order to apply them to any target direction within your field.

.. note::
    In order to run this pipeline you need to have performed `prefactor`_ v3.0 before and provide the corresponding solution set ``solutions.h5``.
    
.. _contents:

Contents
--------

.. toctree::
   :maxdepth: 3
   
   installation
   preparation
   pipelines
   acknowledgments

.. _prefactor: https://github.com/lofar-astron/prefactor
