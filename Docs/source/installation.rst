.. index:: Installation

=================================================
Installation
=================================================

Hardware Requirements
^^^^^^^^^^^^^^^^^^^^^

The pipeline is implemented in the `genericpipeline`_ framework. It is designed to run non-interactively on a cluster via submission to a job queue. The pipeline has been tested on the following type of computing node:

* 2 socket x 16 core (32 threads) 2.10 GHz
* 192 GB RAM
* FDR Infiniband
* 100 TB disk space

For basic pipeline profiling, please see Appendix A in `Morabito`_ et al. (2022). While the configuration can be adapted to your particular cluster specifications, **we recommend at least 32 cores and 192 GB RAM**. Larger number of cores will help reduce the runtime of the pipeline.

The total data volume will reach about 2.5 times that of the raw dataset downloaded from the LTA. If the data is dysco compressed, it will be between 4-6 TB (depending on the number of international stations participating) meaning you will need 10 - 15 TB available. A pre-dysco compression dataset will be around 20 TB and you will need about 50 TB of available disk space. 

.. note::
    Do not forget to check whether your data is dysco compressed! When you stage your data at the LTA you will get a summary of how big it will be.  You will need 2.5 times this size in disk space.

Software Requirements -- with Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the recommended method to run the pipeline. You will need the following the following:

   * An appropriate Singularity image. You may use another one but be aware that there may be software compatibility issues.
     We recommend::
    $ wget https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/lofar_sksp_v3.5_x86-64_generic_noavx512_ddf.sif?action=show

   * The lofar-vlbi github repository (master branch)::

    $ git clone https://github.com/lmorabit/lofar-vlbi

   * The prefactor github repository (see note about aoflagger)::

    $ git clone https://github.com/lofar-astron/prefactor

   * The facet self-cal github repository::

    $ git clone https://github.com/rvweeren/lofar_facet_selfcal

   * The lofar_helpers github repository::

    $ git clone https://github.com/jurjen93/lofar_helpers


   * **NOTE ABOUT AOFLAGGER**
     The latest recommended Singularity version (listed above) has a newer version of aoflagger which does not work with the RFI flagging strategy included in `prefactor`_ and the `LOFAR-VLBI`_ pipeline. You will need to point your parsets to the RFI flagging strategy found here: https://git.astron.nl/eosc/prefactor3-cwl/-/blob/master/rfistrategies/lofar-default.lua . Download this strategy and update the RFI strategy in your prefactor and LOFAR-VLBI parsets to point to your local copy. 



Software Requirements -- without Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If for some reason you are not able to use Singularity, please contact us for instructions. 

.. _genericpipeline: https://www.astron.nl/citt/genericpipeline/
.. _facetselfcal: https://github.com/rvweeren/lofar_facet_selfcal
.. _lofar_helpers: https://github.com/jurjen93/lofar_helpers
.. _Morabito: https://ui.adsabs.harvard.edu/abs/2022A%26A...658A...1M/abstract
.. _Singularity: https://sylabs.io/guides/3.6/user-guide/
.. _LOFAR-VLBI: https://github.com/lmorabit/lofar-vlbi
.. _LoTSS catalogue server: https://vo.astron.nl/lofartier1/lofartier1.xml/cone/form
.. _LBCS catalogue server: https://lofar-surveys.org/lbcs.html
.. _Long Baseline Pipeline GitHub issues: https://github.com/lmorabit/lofar-vlbi/issues
.. _prefactor: https://github.com/lofar-astron/prefactor
.. _prefactor documentation: https://www.astron.nl/citt/prefactor/
.. _prefactor tutorial: https://www.astron.nl/lofarschool2018/Documents/Thursday/prefactor_tutorial.pdf
.. _documentation: file:///media/quasarfix/media/cep3/prefactor/docs/build/html/parset.html
.. _ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline

