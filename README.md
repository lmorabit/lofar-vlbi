Welcome to the LOFAR long baseline pipeline!

This pipeline is designed to be fully compatible with prefactor V3.0.  It uses and produces solutions in h5parm format.

Before running this pipeline, you should have already prefactor on your calibrator and target. The prefactor repository, along with the relevant parsets, can be found here: https://github.com/lofar-astron/prefactor . Please be advised that you should only use 3C 147 or 3C 48 as your standard flux calibrators. If this is not the case, you can copy bandpass solutions from a reference 3C 48 observation (but this is still experimental).

**INSTRUCTIONS FOR PREFACTOR**

1. Run Pre-Facet-Calibrator.parset with the following modification:

Replace

! process_baselines_cal = [CR]S*& 

with

! process_baselines_cal = * 

2. Run Pre-Facet-Target.parset with the default settings.

For both of these steps, you will have to modify things like the path to your data, etc.  Please see the prefactor repository for instructions (https://github.com/lofar-astron/prefactor).

**INSTRUCTIONS FOR THE LONG BASELINE PIPELINE**

1. Run LB-Delay-Calibrator.parset

Please update the necessary parameters in the "Please update these parameters" section of the parset.

*This will work out of the box if your field is already covered by both LBCS and LoTSS. If this is not the case, the aumatic catalogue generation step will fail. This can be fixed by providing (a) manual catalogue(s) with the right format.*

After this step, the data will have all the prefactor solutions applied, in the DATA column. The CORRECTED_DATA column is DATA + phase and amplitude solutions from the best in-field calibrator. However, the next step starts with the DATA column. 


2. Run LB-Split-Calibrators.parset

Please update the necessary parameters in the "Please update these parameters" section of the parset.

*If you have used a non-standard catalogue, please either name it the same as in the "These parameters may need to be updated" section of the parset, or change the name of the delay_cat there.*

After this step, you will have a number of smaller measurement sets with phased-up core stations combined into ST001, for all the directions of LBCS calibrators.  These will be accompanied by the TEC solutions in each direction. 1

