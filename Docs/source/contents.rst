.. lofar-vlbi documentation master file, created by
   sphinx-quickstart on Mon May 11 11:37:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for the LOFAR-VLBI pipeline
=========================================

The LOFAR-VLBI pipeline is a calibration and imaging pipeline which includes all of LOFAR's international stations. It looks and solves for delays using an in-field calibrator within the field-of-view of the pointing. The pipeline provides solution tables (in h5parm format) and a self-calibrated, corrected dataset for this in-field calibrator. The delay solutions are applied back to the original data, and from there the pipeline can split out smaller datasets in the desired target direction(s) within the field of view.

Table of Contents
^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   installation
   license
   acknowledgements

Quick Links
^^^^^^^^^^^

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
