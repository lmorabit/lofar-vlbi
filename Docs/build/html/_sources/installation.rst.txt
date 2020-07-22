.. index:: Installation

=================================================
Installation
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

Software Requirements -- with Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the recommended method to run the pipeline. You will need:
* The recommended Singularity image
* The lofar-vlbi github repository (master branch)
* The prefactor github repository (commit 7e9103d10c8e37ee2ac2203b678af295ed03e4fd)

and that's it! Specific instructions on how to get these are below.

**Singularity**

The LOFAR-VLBI pipeline and prefactor (specific commit listed above) have been successfully tested against a specific singularity image. You can download this image by running::

       $ wget https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/lofar_sksp@e7512b7e92e5a3124e6ed177650e99a8c4eb2263_with_pyvo.sif

The singularity image was built with singularity version 2.5.  Users have reported that it also works with singularity versions 3.1 and 3.3.

**LOFAR-VLBI pipeline**

To get the LOFAR-VLBI pipeline, run::

    $ git clone https://github.com/lmorabit/lofar-vlbi

**Prefactor**

To get the specific version of prefactor you will need, run::

    $ git clone https://github.com/lofar-astron/prefactor
    $ cd prefactor
    $ git checkout 7e9103d10c8e37ee2ac2203b678af295ed03e4fd

Software Requirements -- without Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If for some reason you are not able to use Singularity, the following software must be installed locally, *in addition to* the LOFAR-VLBI and Prefactor github repositories listed above:

   * AOFlagger >= 2.14.0
   * DP3 >= 4.0
   * LoSoTo >= 2.0
   * Montage 
   * prefactor (commit 7e9103d10c8e37ee2ac2203b678af295ed03e4fd)
   * PGPLOT
   * Difmap (built with modified version of corplt.c found here: https://github.com/nealjackson/loop3_difmap/corplt.c)

Along with the following python packages:
   * packages listed in requirements.txt from https://github.com/lmorabit/lofar-vlbi
   * PyBDSF
   * RMextract

The PyBDSF and RMextract python packages can be installed either by cloning the github repositories, or via pip::

   $ pip install https://github.com/lofar-astron/PyBDSF/archive/v1.9.1.tar.gz
   
   $ git clone https://github.com/lofar-astron/RMextract
   $ pip install -e RMextract

The other packages can be installed by running::

   $ pip install -r requirements.txt

The following software is **optional**, but must be installed for their respective features to work:

   * The ddf-pipeline (https://github.com/mhardcastle/ddf-pipeline) and its software prerequisites (listed in docs/manual.md)

.. _genericpipeline: https://www.astron.nl/citt/genericpipeline/
