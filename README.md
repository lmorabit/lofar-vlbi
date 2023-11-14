

** **NOTE: The genericpipeline version is no longer supported but can be found as tag v5.0-genericpipeline and branch v5.0-genericpipeline-branch** **

[![Documentation Status](https://readthedocs.org/projects/lofar-vlbi/badge/?version=latest)](https://lofar-vlbi.readthedocs.io/en/latest/?badge=latest)

More complete documentation can be found on https://lofar-vlbi.readthedocs.io/en/latest/


This is the LOFAR VLBI data reduction repository. It contains software and workflows that can be used to process LOFAR data to high (sub-arcsecond to arcsecond) resolution. 

It has dependencies on the following repositories:
- https://github.com/rvweeren/lofar_facet_selfcal
- https://github.com/jurjen93/lofar_helpers
- https://github.com/tikk3r/flocs/
- https://git.astron.nl/RD/LINC

**Getting started**

The pipeline is built in a modular way, where the desired outcome can be achieved by stringing the different modules together. The modules should be interoperable, but this is still in development to ensure the "glue" between modules works every time.

__Module 1: Pipeline preparation__

***Required** before any pipeline run.*
- Script: plot_field.py 
- Dependencies: VLASS_dyn_summary.php
- Requirements: pyvo, casacore.tables, astropy, numpy, time, requests

What it does: Queries online databases for LotSS and LBCS to construct catalogues of sources in the field, and provide an initial graphical representation. The catalogues that this step generates are **required** for the pipeline to run, so this step cannot be skipped. It will also provide a short summary of the observation parameters, to help the user understand if the data is suitable for LOFAR-VLBI processing. *This step can and should be run before downloading the data to avoid wasting time and resources on an unsuitable dataset!*

__Module 2: LINC/Prefactor__

***Required** before any pipeline run.*

*Option 1:* Suitable prefactor3 or LINC solutions exist for the target field, i.e. the international stations were included in the calibrator run. In this case, all you need is the final cal_solutions.h5 from the target run (which contains the calibrator solutions). Check that the solutions are good before proceeding. 

*Option 2:* Suitable calibrator / target solutions do not already exist. You need to run LINC (https://git.astron.nl/RD/LINC) for the calibrator and the target to generate the final cal_solutions.h5 from the target run (which contains the calibrator solutions). Do not change any of the default settings for either the calibrator or target, but please do check that your solutions are good before proceeding.













