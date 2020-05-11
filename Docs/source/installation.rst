======================================
Installation and software requirements
======================================

Software Requirements
^^^^^^^^^^^^^^^^^^^^^

To get the LOFAR-VLBI pipeline, run::

    $ git clone https://github.com/lmorabit/lofar-vlbi

You will also need a local copy of prefactor, run::

    $ git clone https://github.com/lofar-astron/prefactor
    $ cd prefactor
    $ git checkout 7e9103d10c8e37ee2ac2203b678af295ed03e4fd

The software dependencies for the LOFAR-VLBI pipeline and prefactor are listed below. 

.. note::
    Fortunately, everything has been packaged into a singularity image, which can be found here:

    https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/lofar_sksp@e7512b7e92e5a3124e6ed177650e99a8c4eb2263.sif?action=show

    This version has been tested, but does not work with the -H flag in the h5exp_gsm step in the prefactor Pre-Facet-Target.parset, so you have to remove that manually. 

If you do not wish to use this singularity image, the following software must be installed locally:

   * AOFlagger >= 2.14.0
   * DP3 >= 4.0
   * LoSoTo >= 2.0
   * Montage 
   * prefactor (commit 7e9103d10c8e37ee2ac2203b678af295ed03e4fd)
   * PGPLOT
   * Difmap (built with modified version of corplt.c found here: https://github.com/nealjackson/loop3_difmap/corplt.c)
   
The following software is **optional**, but must be installed for their respective features to work:

   * The ddf-pipeline (https://github.com/mhardcastle/ddf-pipeline) and its software prerequisites (listed in docs/manual.md)


Python packages
^^^^^^^^^^^^^^^

These are all included in the singularity image, which we strongly recommend using. If you wish to run the pipeline locally, the required Python packages are listed in requirements.txt. One can easily install these using::

   pip install -r requirements.txt

The exceptions are:

   * PyBDSF
   * RMextract

These can both be installed from their respective repositories, or with pip via::

   pip install https://github.com/lofar-astron/PyBDSF/archive/v1.9.1.tar.gz
   
   git clone https://github.com/lofar-astron/RMextract
   pip install -e RMextract
