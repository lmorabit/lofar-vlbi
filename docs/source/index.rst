.. long_baseline_pipeline documentation master file, created by
   sphinx-quickstart on Mon Nov 11 17:22:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

long_baseline_pipeline: a calibration pipeline for LOFAR's international stations
=================================================================================
The LOFAR **lofar-vlbi** analyses LOFAR data incorporating all international stations. It looks and solves for delays using an in-field calibrator within the field-of-view of the pointing. The pipeline provides solution tables (in h5parm format) and a self-calibrated, corrected dataset for this in-field calibrator. The delay solutions are applied back to the original data, and from there the pipeline can split out smaller datasets in the desired target direction(s) within the field of view. 

.. note::
    Before running this pipeline, you need prefactor solutions and, optionally, ddf-pipeline solutions.  Please see the preparation section for more details. 
.. note::
   We highly recommend the use of singularity to run this pipeline (and prefactor). More details are in the preparation section. 
    
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
