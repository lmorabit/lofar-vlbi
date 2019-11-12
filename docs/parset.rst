#################################
Direction independent calibration
#################################
This section will describe the steps performed by ``LB-Delay-Calibration.parset``. Within the parset, the following major steps are performed:

   1. Initial setup
   2. Target preparation
   3. A-team clipping
   4. AOFlagging [toggle]
   5. Subband concatenation
   6. Application of ddf-pipeline DI solutions  [toggle]
   7. Self calibration on the best infield calibrator   [toggle]
   8. Application of self calibration solutions [toggle]

Options marked with ``[toggle]`` can be turned on of off by the user. Other steps are always exectuted.

=============
Initial setup
=============
Initial setup consists of the following steps::

   mk_results_dir
   mk_inspect_dir
   mk_cal_values_dir
   createmap_target
   createmap_target_list
   cp_cal_solutions
   download_cats

``mk_results_dir``: creates the directory where results will be stored.

``mk_inspect_dir``: creates the directory where inspection plots will be stored.

``mk_cal_values_dir``: creates the directory where calibration tables will be stored.

``createmap_target``: creates a mapfile for the target subbands.

``createmap_target_list``: creates a mapfile containing all target subbands as one big list.

``cp_cal_solutions``: copies the prefactor calibration solutions to the inspection directory.

``download_cats``: downloads LoTSS and LBCS catalogues for the field.
