.. index:: Installation

=================================================
Installation, software, and hardware requirements
=================================================

Hardware Requirements
^^^^^^^^^^^^^^^^^^^^^

The pipeline is implemented in the `genericpipeline`_ framework. It is designed to run non-interactively on a cluster via submission to a job queue. Some of the steps in the pipeline require interaction between all running processes, and currently this means the pipeline is limited to running on a single node. The pipeline has been tested on the following type of computing node:

* 2 socket x 16 core (32 threads) 2.10 GHz
* 192 GB RAM
* FDR Infiniband
* 100 TB disk space

With these specifications, the two steps of the pipeline (Delay-Calibration and Split-Directions) will take about 7-8 days and 1-3 days each. The first step can be shortened to about 2-3 days if aoflagging and A-team clipping is turned off. While the configuration can be adapted to your particular cluster specifications, **we recommend at least 32 cores and 192 GB RAM**. Larger number of cores will help reduce the runtime of the pipeline.

The total data volume will reach about 2.5 times that of the raw dataset downloaded from the LTA. If the data is dysco compressed, it will be between 4-6 TB (depending on the number of international stations participating) meaning you will need 10 - 15 TB available. A pre-dysco compression dataset will be around 20 TB and you will need about 50 TB of available disk space. 

.. note::
    Do not forget to check whether your data is dysco compressed! When you stage your data at the LTA you will get a summary of how big it will be.  You will need 2.5 times this size in disk space.

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

    https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/lofar_sksp@e7512b7e92e5a3124e6ed177650e99a8c4eb2263_with_pyvo.sif

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


.. _genericpipeline: https://www.astron.nl/citt/genericpipeline/
