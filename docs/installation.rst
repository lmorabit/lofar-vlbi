################################################
Obtaining and running the long baseline pipeline
################################################
.. contents::
      :local:

**********************
Obtaining the pipeline
**********************

To obtain the latest version of the long baseline pipeline, run::

   git clone https://github.com/lmorabit/long_baseline_pipeline.git

This will pull the latest master branch.

.. contents::
      :local:

Software requirements
=====================
The following software _must_ be installed, in order to run the pipeline:

   * AOFlagger
   * DP3 >= 4.0
   * LOFAR software
   * LoSoTo >= 2.0
   * prefactor >= 3.0

The following software is _optional_, but must be installed for their respective features to work:

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

Environment settings
====================
The generic pipeline framework needs to be able to find the scripts. This is set in the config file that is passed to the generic pipeline. In order to achieve this, the path to the long baseline pipeline needs to be added to the recipe_directories list in this file.

********************
Running the pipeline
********************

.. contents::
      :local:

Once all parameters are set, the pipeline can be run as, for example,::

   genericpipeline.py -c pipeline.cfg LB-Delay-Calibration.parset

Below setup of the parset is outlined.

Configuring the LB-Delay-Calibration parset
===========================================
This sections describes how to set up the ``LB-Delay-Calibration`` parset. Only sections that are relvant for the user to change are discussed.

Software
--------
::

   # software paths
   ! long_baseline_pipeline_dir   = /path/to/long_baseline_pipeline
   ! prefactor_dir                = /path/to/prefactor
   ! losoto_directory             = /path/to/losoto
   ! aoflagger		          = /path/to/aoflagger/bin/aoflagger ## path to your aoflagger executable
   ! lofar_directory              = $LOFARROOT

This section points the pipeline to the required software. Each should point to the respective software's base directory, except for aoflagger, which should point to the executable. ``lofar_directory`` should point to the base directory of the lofar software, which is what ``$LOFARROOT`` typically points to in a standard environment. Hence this usually will not need changing.

Observation data
----------------
::

   # ## target data information
   ! target_input_path    = /data/scratch/lb_bw/targetdata
   ! target_input_pattern = L*.MS

The pipeline will look for the raw data under ``target_input_path``, where it looks for Measurement Sets that follow the given ``target_input_pattern``.

Prefactor solutions
-------------------
::

   ## Prefactor solution information
   ! cal_solutions = /path/to/prefactor/solutions.h5
   ! solutions	   = input.output.job_directory/solutions.h5
   ! cal_table	   = combinedsols
   ! phasesol      = TGSSphase

Solutions from the prefactor pipeline are pre-applied before further calibrating the international stations. ``cal_solutions`` points to the final H5parm produced after running prefactor's target pipeline. ``cal_table`` is the name of the output solset after adding the international stations to the prefactor target phase solutions, which are stored in the ``phasesol`` soltab.

Averaging and flagging
----------------------
::

   ## Stations to flag
   ! flag_baselines         = [ ] ## for HBA data before October 2015, should set to: [ CS013HBA* ]

   ## averaging information -- do not touch unless you know what you are doing!
   ! cal_shift_avg_freqstep = 4
   ! cal_shift_avg_timestep = 4

Here the user can set the averaging parameters for LBCS calibrators that are split off. Using ``flag_baselines``, the user can explicitely specify baselines and/or stations that need to be flagged. For the sytax, see the DP3 documentation. ``cal_shift_avg_freqstep`` is the factor with which to average in frequency. ``cal_shift_avg_timestep`` is the factor with which to average in time.

DDF pipeline options
--------------------
::

   ! substep_ddf  = # ,apply_ddf   ## set to apply_ddf to use, leave blank otherwise
   ! ddf_soldir	  = /data/scratch/lb_bw/targetddf/L602888/SOLSDIR	 ## path to 'SOLSDIR' from ddf-pipeline run
   ! delaycal_col = DATA ## set to DATA_DI_CORRECTED if applying ddf solutions, DATA otherwise.

If the field has been processed by the ddf-pipeline, it is possible to apply the DIS2 solutions, in order to subtract the 6'' map later on. This can be enabled by setting ``substep_ddf`` to ``,apply_ddf`` (the comma is important). ``ddf_soldir`` then points to the SOLSDIR of the ddf-pipeline run of the field. Finally, ``delaycal_col`` should be set to either DATA_DI_CORRECTED, when applying these solutions, or DATA if not using ddf-pipeline solutions.

Substeps
--------
::

   ## setup for running steps
   ! substep_clip  = clipATeam, 	## leave blank for no a-team clipping
   ! substep_flag  = # aoflagging,	## aoflagging,  ## turned off by default
   ! substep_delay = # ,delaycal	## off by default, turn on to do a single in-field delay calibration
   ! substep_apply = # ,applydelay ## off by default, turn on ONLY IF you have run substep_delay and wish to apply solutions

These options control various substeps that the pipeline will do. In total, there are four steps that can be enabled: A-team clipping, AOFlagging, delay calibration and applying delay calibration solutions. Each step can be enabled by removing the leading # (but keeping the comma), or disabled leaving it empty. Enabling ``substep_delay`` will search for the best infield calibrator, followed by self calibration on it. If ``substep_apply`` is enabled, these solutions are applied to the data. The calibrated data is then stored in the CORRECTED_DATA column.
