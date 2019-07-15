Welcome to the LOFAR long baseline pipeline!

This pipeline is designed to be fully compatible with prefactor V3.0.  It uses and produces solutions in h5parm format.

Before running this pipeline, you should have already prefactor on your calibrator and target. The prefactor repository, along with the relevant parsets, can be found here: https://github.com/lofar-astron/prefactor . Please be advised that you should only use 3C 147 or 3C 48 as your standard flux calibrators. If this is not the case, you can copy bandpass solutions from a reference 3C 48 observation (but this is still experimental, and not yet built into this pipeline). 

**Software requirements**
* Prefactor: https://github.com/lofar-astron/prefactor
* LoSoTo: https://github.com/revoltek/losoto
* aoflagger: https://sourceforge.net/p/aoflagger/wiki/Home/
* **_Optional_**: ddf-pipeline: https://github.com/mhardcastle/ddf-pipeline

**Getting help**
* If you have issues with the installation or running of the software requirements, please address these externally from this repository, e.g., if you have a problem with prefactor please open an issue on the prefactor repository.
* If you have trouble running any of the parsets from this repository, please _open an issue_ rather than emailing! If you are going to submit an issue, please check those already open -- someone else may have already had the same problem.
* For members of the Long Baseline Working Group, there is a private slack channel. Please email the working group for details.

**INSTRUCTIONS FOR PREFACTOR**

1. Run Pre-Facet-Calibrator.parset with the following modification:

Replace  
> ! process_baselines_cal = [CR]S*&   

with  

> ! process_baselines_cal = * 

2. Run Pre-Facet-Target.parset with the default settings.

For both of these steps, you will have to modify things like the path to your data, etc.  Please see the prefactor repository for instructions (https://github.com/lofar-astron/prefactor).

**INSTRUCTIONS FOR DDF**

*This is an **optional** step which is only relevant if you wish to progress to wide-field imaging afterwards, and is not necessary to run the pipeline.*

1. Run the ddf-pipeline (following instructions on their github repository) with the LOFAR Surveys default settings. This will operate on the output of Pre-Facet-Target and only provides additional phase solutions for the core and remote stations. 

**INSTRUCTIONS FOR THE LONG BASELINE PIPELINE**

___Notes:___ The pipeline solves for TEC, which is a frequency dependent effect. You need to process an absolute minimum of 10 MHz of bandwidth (30 subbands is safe) for the pipeline to run. 

1. Run LB-Delay-Calibrator.parset  

Please update the necessary parameters in the "Please update these parameters" section of the parset.  
_optional_: If you have run the ddf-pipeline, please update the DDF options section as well.

*This will work out of the box if your field is already covered by both LBCS and LoTSS. If this is not the case, the aumatic catalogue generation step will fail. This can be fixed by providing (a) manual catalogue(s) with the right format.*

After this step, the data will have all the prefactor solutions applied, in the DATA column. If you have applied ddf-pipeline solutions, they will be in the CORRECTED_DATA column. If you optionally choose to run the delay calibration and apply steps (not necessary for the next step) then you will also have new measurement sets with those corrections applied.

2. Run LB-Split-Calibrators.parset 

Please update the necessary parameters in the "Please update these parameters" section of the parset.  
_optional_: If you have run the ddf-pipeline, please update the DDF options section as well.

*If you have used a non-standard catalogue, please either name it the same as in the "These parameters may need to be updated" section of the parset, or change the name of the delay_cat there.*

After this step, you will have a number of smaller measurement sets with phased-up core stations combined into ST001, for all the directions of LBCS calibrators.  These will be accompanied by the TEC solutions in each direction, and loop3 will have been run on new measurement sets where TEC has been applied.

