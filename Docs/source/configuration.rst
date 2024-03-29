.. index:: Configuration

*************
Configuration
*************

===================
Running singularity
===================

Depending on the installation of `Singularity`_ you have access to, the syntax may change slightly. The following instructions are for Singularity version 3.1. 

You can either start a Singularity shell, which you can then use interactively just as you would any other shell (e.g., bash, csh, etc.), or you can execute a script using Singularity. In both cases, you need to:

 * provide an image for Singularity to use
 * 'bind' directories that will need to be used

For example, to start a Singularity shell::

        $ singularity shell --bind /scratch/data/,/home/username/software my_singularity_image.sif

This will return a prompt where you will be 'inside' the Singularity image -- i.e., have access to everything in it.

To run a script (which is generally how you will run the pipeline), you would do::

        $ singularity exec --bind /scratch/data,/home/username/software my_singularity_image.sif myscript.py arg1

These are just examples; please refer to the documentation for your specific version of `Singularity`_ for further information.

=================================
Running a Generic Pipeline parset
=================================

Both `prefactor`_ and the `LOFAR-VLBI`_ pipeline use the `genericpipeline`_ framework. The instructions for a pipeline are written in a specific format in a *parset*, which is passed along with a configuration file, to the generic pipeline. For example::

   $ genericpipeline.py -c pipeline.cfg pipeline.parset

Each parset will contain different instructions, relating to what that part of the pipeline is doing. The `prefactor`_ parsets consist of:

 * Pre-Facet-Calibrator.parset
 * Pre-Facet-Target.parset

And the `LOFAR-VLBI`_ parsets are:

 * Delay-Calibration.parset
 * Split-Directions.parset (still under development)

These can all use the same configuration file (see next section).

===========================
Pipeline configuration file
===========================

If you are using the recommended Singularity image, the configuration file is relatively straightforward to set up. The *pipeline.cfg* file from the `LOFAR-VLBI`_ github repository is already configured for the Singularity image. You only need to change the following lines::

        runtime_directory=
        recipe_directories=[%(pythonpath)s/lofarpipe/recipes,/opt/lofar/losoto] 

Set the ``runtime_directory`` to wherever your parent directory is for the pipeline to run; e.g if you set this to ``runtime_directory=/scratch/myusername/`` then running the *Pre-Facet-Calibrator.parset* will create the directory ``/scratch/myusername/Pre-Facet-Calibrator`` and run in that directory.

To the ``recipe_directories`` list you need to add where you cloned the `prefactor`_ and `LOFAR-VLBI`_ repositories, e.g. if you put them in ``/scratch/`` then::

        recipe_directories=[%(pythonpath)s/lofarpipe/recipes,/opt/lofar/losoto,/scratch/prefactor,/scratch/lofar-vlbi] 

The settings at the bottom of the configuration file::

        [remote]
        method = local
        max_per_node = 32

Are sensible for running on a single node with 32 cores available. Please consult your system administrator if you think this needs to be changed.  

If you are using a local version of the LOFAR software instead of the Singularity image, please consult your system administrator for help on setting up the configuration file.

========================
Using your own catalogue
========================

The `LOFAR-VLBI`_ pipeline will automatically try to download information from both the `LBCS catalogue server`_ and the `LoTSS catalogue server`_. Both of these are required to help select the best in-field calibrator. As described in the `Before you begin`_ section, if LoTSS coverage for your field does not yet exist, you can manually make your own field catalogue, or if you have information on calibrator sources from somewhere else, you can provide your own catalogue.  This needs to be a CSV file with a minimum of columns **Source_name,RA,DEC** with your in-field calibrator as the top entry. The columns must be named as described here; you can have more columns than this but **Source_name,RA,DEC** must exist. 

If you use your own catalogue, update the following line in **Delay-Calibration.parset**::

    ! lotss_skymodel         = {{ results_directory }}/lotss_catalogue.csv

to the absolute path for your csv file. It does not need to be named lotss_catalogue.csv.  You do not need to make any further changes to the catalogue.

If there is no LBCS coverage for your field, please contact someone from the LOFAR-VLBI working group.

===============================
Setting the directions to image
===============================

The **Delay-Calibration** step generates some output catalogues, which are stored in its *results* directory. These include:

* delay_calibrators.csv - a list of potential LBCS calibrators in the field 
* image_catalogue.csv - everything else

Once the **Delay-Calibration** step has run, you can simply edit or replace the *image_catalogue.csv* file to include only the source(s) you wish to image. The more directions you want to image, the longer the pipeline will take, so you should really limit this to your target of interest. The file needs to be in **csv format** with the **same column names** as *image_catalogue.csv* and flux densities in Janskys.

Selecting imaging parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The imaging parameters for the delay calibrator are stored in the *facetselfcal_config.txt* file. These should not be changed unless you know what you are doing. 
   
.. _help:


.. _Before you begin: before.html
.. _genericpipeline: https://www.astron.nl/citt/genericpipeline/
.. _Singularity: https://sylabs.io/guides/3.6/user-guide/
.. _LOFAR-VLBI: https://github.com/lmorabit/lofar-vlbi
.. _LoTSS catalogue server: https://vo.astron.nl/lofartier1/lofartier1.xml/cone/form
.. _LBCS catalogue server: https://lofar-surveys.org/lbcs.html
.. _Long Baseline Pipeline GitHub issues: https://github.com/lmorabit/lofar-vlbi/issues
.. _prefactor: https://github.com/lofar-astron/prefactor
.. _prefactor documentation: https://www.astron.nl/citt/prefactor/
.. _documentation: file:///media/quasarfix/media/cep3/prefactor/docs/build/html/parset.html
.. _ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
