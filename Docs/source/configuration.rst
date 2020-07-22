.. index:: Configuration

*************
Configuration
*************

=============================
Preparing the data: Prefactor
=============================

The LOFAR-VLBI pipeline makes use of prefactor solutions to apply to the data. Therefore you must pre-process your data through prefactor, both the calibrator and target pipelines. For instructions how to run prefactor on your data, please look at the `prefactor documentation`_. For any issues you encounter with prefactor, please open an issue on the `prefactor`_ repository.


* The default is now for **Pre-Facet-Calibrator.parset** to process all stations, including the international stations. You can run this with the default settings. Please check the outputs to make sure they are sensible! 

.. note::
    If your standard calibrator is either 3C 295 or 3C 196, the standard models in the prefactor github repository do not have sufficiently high resolution, but high-resolution models do exist. Please contact the long baseline working group for help. 

* The **Pre-Facet-Target.parset** should be run with all the standard defaults. This will copy over the solutions from Pre-Facet-Calibrator and add the self-cal phase solutions for the core and remote stations, which are necessary for the LOFAR-VLBI pipeline. Please check the outputs to make sure they are sensible!  Also note any stations which were flagged as 'bad' as you will need to pre-flag these for the LOFAR-VLBI pipeline.

.. note::
    Processing of interleaved datasets is not currently supported.

======================
Optional: ddf-pipeline
======================

This is an optional step and is not necessary to run the pipeline unless you do not have a LoTSS-style catalogue for your field. The ddf-pipeline requires some advanced user knowledge to set up and run, so please contact lofar-admin@strw.leidenuniv.nl if you are considering doing this step. If you are using Surveys data it may have already been run for your pointing; if not, we can help. 
Collaborative projects with the Surveys KSP are also possible, if you have your own data and want it processed through the SKSP infastructure to carry out this step. Contact lofar-admin@strw.leidenuniv.nl for more details. 


This step is only necessary in the case where your field has not been covered yet by LoTSS, to generate a catalogue of sources in the field which is used by the LOFAR-VLBI pipeline to help select the best candidate for in-field calibration. If you can query sources in your field with the `LoTSS catalogue server`_ then you do not need to generate this catalogue. 

.. note::
    The recommended singularity image works with prefactor and the LOFAR-VLBI pipeline, but not the ddf-pipeline.  Please refer to the `ddf-pipeline`_ documentation for its separate software requirements, or contact lofar-admin@strw.leidenuniv.nl .


The `ddf-pipeline`_  operates on the results of Pre-Facet-Target and provides:

* additional phase solutions for core and remote stations
* a self-calibrated image at 6" resolution
* an initial catalogue of sources in the field

To generate the final catalogue, use the *quality_pipeline.py* script (found in the `ddf-pipeline`_ *scripts* sub-directory) with an appropriate configuration file (the example *quality-example.cfg* is in the *examples* sub-directory). The bootstrap catalogues can be downloaded from here: https://www.extragalactic.info/bootstrap/ . Note that you will also need to convert all your flux values from Jy to mJy.

The LOFAR-VLBI pipeline **requires** the information on the sources, either from this output catalogue or the `LoTSS catalogue server`_ , and if you run the ddf-pipeline it can use the additional phase solutions (but this is not required). We recommend skipping this step if your field is already in the `LoTSS catalogue server`_ unless you wish do do wide-field imaging at high resolution, rather than imaging science targets in a few (or one) directions. 


===============================
Running the LOFAR-VLBI pipeline
===============================

The LOFAR-VLBI pipeline uses the same ``genericpipeline`` framework as prefactor. You can see the prefactor `documentation`_ on how to modify the ``pipeline.cfg`` and the corresponding parset files before you start the pipeline, although you should already be familiar with this if you've done it for prefactor.

.. note::
    The pipeline.cfg file in the `LOFAR-VLBI`_ repository already contains paths for the singularity image, although some paths will need to be local. Please check this file carefully before making changes. 

The LOFAR-VLBI pipeline is broken into two steps: **Delay-Calibration.parset** and **Split-Directions.parset**. The first parset does all the heavy lifting; it applies the prefactor solutions, splits out best in-field calibrator candidate, performs the delay calibration on it, and applies these corrections back to the data. The second parset takes the resulting CORRECTED_DATA, splits out the directions in which you wish to image, and runs self-calibration on them. 


Before running the pipeline, you should check:

* If there are any bad stations flagged by prefactor. These will need to be manually input into the parsets. Follow exactly the syntax for the example given in the parset.

* Check the rest of the "Please update these parameters" section. Comments in the parset(s) describe what they are. 

* Optional: if you have run the ddf-pipeline, please update the DDF options as well. If you are only using the catalogue, update the lotss_skymodel parameter to point to your output file. 

Once all parameters are set, the pipeline can be run as, for example::

   genericpipeline.py -c pipeline.cfg Delay-Calibration.parset

========================
Using your own catalogue
========================

The pipeline will automatically try to download information from both the `LBCS catalogue server`_ and the `LoTSS catalogue server`_. Both of these are required to help select the best in-field calibrator. You can generate an appropriate catalogue to replace the LoTSS catalogue by running the `ddf-pipeline`_ and then the *quality_pipeline.py* script. The output catalogue will be named *image_full_ampphase_di_m.NS.cat.fits*.  The only thing you need to do is convert this to a csv file, and then update the following line in **Delay-Calibration.parset**::

    ! lotss_skymodel         = {{ results_directory }}/lotss_catalogue.csv

to the absolute path for your csv file. It does not need to be named lotss_catalogue.csv.  You do not need to make any further changes to the catalogue.

If there is no LBCS coverage for your field, please contact someone from the LOFAR-VLBI working group.

===============================
Setting the directions to image
===============================

The **Delay-Calibration** step generates some output catalogues, which are stored in its *results* directory. These include:

* delay_calibrators.csv - a list of potential LBCS calibrators in the field 
* best_delay_calibrators.csv - the best LBCS calibrator to use for the delay calibration
* subtract_sources.csv - bright sources and LBCS calibrators that may need to be subtracted to improve image fidelity
* image_catalogue.csv - everything else

Once the **Delay-Calibration** step has run, you can simply edit or replace the *image_catalogue.csv* file to include only the source(s) you wish to image. The more directions you want to image, the longer the pipeline will take, so you should really limit this to your target of interest. The file needs to be in **csv format** with the **same column names** as *image_catalogue.csv* and flux densities in Janskys.

Selecting imaging parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the pipeline will run self-calibration using difmap. This is an order of magnitude faster (usually ~30 min) than any self-calibration using native LOFAR tools, and already optimised for VLBI. Difmap operates on the XX and YY polarisations independently, but the self-calibration script converts these solutions to an h5parm, applies them, and makes a Stokes I image from the corrected data using wsclean. The final self-calibrated dataset will have TEC-corrected, un-self-calibrated data in the **DATA** column and TEC + self-cal corrected data in the **CORRECTED_DATA** column. The user is free to perform more self-calibration, or re-do the self-calibration, using any tools they wish. The data at this point is already corrected for beam effects (including the array factor), so you are free to use any imaging / gain calibration software you like.

The self-calibration script run by the pipeline has the following default parameters:
* Number of pixels = 512
* Pixel scale = 50 milli-arcsec

This gives an image which is 25.6 x 25.6 arcseconds. If your source is larger than this, you will need to adjust the number of pixels, following the convention of using powers of 2 (512,1024,2048,... etc.). 
   
.. _help:

.. _LOFAR-VLBI: https://github.com/lmorabit/lofar-vlbi
.. _LoTSS catalogue server: https://vo.astron.nl/lofartier1/lofartier1.xml/cone/form
.. _LBCS catalogue server: https://lofar-surveys.org/lbcs.html
.. _Long Baseline Pipeline GitHub issues: https://github.com/lmorabit/lofar-vlbi/issues
.. _prefactor: https://github.com/lofar-astron/prefactor
.. _prefactor documentation: https://www.astron.nl/citt/prefactor/
.. _documentation: file:///media/quasarfix/media/cep3/prefactor/docs/build/html/parset.html
.. _ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline
