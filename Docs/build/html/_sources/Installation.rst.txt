================================================
Obtaining and running the long baseline pipeline
================================================


Obtaining the pipeline
======================

To obtain the latest version of the long baseline pipeline, run::

   git clone https://github.com/lmorabit/long_baseline_pipeline.git

This will pull the latest master branch.

Software requirements
=====================
The following software **must** be installed, in order to run the pipeline:

   * AOFlagger
   * DP3 >= 4.0
   * LOFAR software
   * LoSoTo >= 2.0
   * Montage
   * prefactor >= 3.0

The following software is **optional**, but must be installed for their respective features to work:

   * DDFacet

Python packages
===============
Most required Python packages are put in requirements.txt. One can easily install these using::

   pip install -r requirements.txt

The exceptions are:

   * python-casacore
   * RMextract

These can both be installed from their respective repositories, or with pip via::

   pip install https://github.com/lofar-astron/PyBDSF/archive/v1.9.1.tar.gz
   
   git clone https://github.com/lofar-astron/RMextract
   pip install -e RMextract
