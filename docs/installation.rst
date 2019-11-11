****
Obtaining the long baseline pipeline
****

To obtain the latest version of the long baseline pipeline, run::

   git clone https://github.com/lmorabit/long_baseline_pipeline.git

This will pull the latest master branch.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Software requirements
#####
The following software _must_ be installed, in order to run the pipeline:

   * AOFlagger
   * DP3 >= 4.0
   * LOFAR software
   * LoSoTo >= 2.0
   * prefactor >= 3.0

The following software is _optional_, but must be installed for their respective features to work:

   * DDFacet

Python packages
####
Most required Python packages are put in requirements.txt. One can easily install these using::

   pip install -r requirements.txt

The exception is:

   * python-casacore
