

** **NOTE: The genericpipeline version is no longer supported but can be found as tag v5.0-genericpipeline and branch v5.0-genericpipeline-branch** **

[![Documentation Status](https://readthedocs.org/projects/lofar-vlbi/badge/?version=latest)](https://lofar-vlbi.readthedocs.io/en/latest/?badge=latest)

More complete documentation can be found on https://lofar-vlbi.readthedocs.io/en/latest/ (UNDER CONSTRUCTION)


This is the LOFAR VLBI data reduction repository. It contains software and workflows that can be used to process LOFAR data to high (sub-arcsecond to arcsecond) resolution. 

It has dependencies on the following repositories:
- https://github.com/rvweeren/lofar_facet_selfcal
- https://github.com/jurjen93/lofar_helpers
- https://github.com/tikk3r/flocs/
- https://git.astron.nl/RD/LINC
- https://git.astron.nl/RD/VLBI-cwl/
- https://github.com/mhardcastle/ddf-pipeline (optional)

**Getting started**

The pipeline is built in a modular way, where the desired outcome can be achieved by stringing the different modules together. The modules should be interoperable, but this is still in development to ensure the "glue" between modules works every time.

The setup and running of LINC and VLBI-cwl modules can be accomplished using the runners in https://github.com/tikk3r/flocs/runners although the submission scripts may need adjustments per cluster. 

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

__Module 3: ddf-pipeline__

***Required** if you are planning on doing widefield processing.*

***Not required** if you are only processing individual postage stamps.*

What it does: This is the direction-dependent processing for the Dutch array. Solutions from this are used to subtract sources outside the field of view of the international stations when performing widefield, high resolution imaging.  If you are only processing "postage stamp" images of individual targets in the field, this step is not required.

__Module 4: VLBI-cwl__

***Required** if you are planning on doing widefield processing.*

This module applies all relevant solutions from the previous modules, carries out the initial in-field calibration on a delay calibrator, and splits out individual target directions for further self-calibration. There are two main workflows:

- delay-calibration.cwl
- split-directions.cwl

The first of these workflows prepares the data and the end result is the first in-field calibration on a calibrator source from the Long Baseline Calibrator Survey (LBCS).  It runs several sub-workflows:

- setup.cwl
- concatenate-flag.cwl
- phaseup-concat.cwl

These can be run directly, with some book-keeping to gather the results of one sub-workflow to input to the next sub-workflow. This is advantageous when submitting jobs to a queuing system with limits on the walltime for jobs, as the entire workflow takes about a week to run. The setup step is the longest. 

The second workflow splits out directions of interest for further self-calibration. 

__Module 5a: full FoV, intermediate resolution (1-2 arcsec) imaging__

PIPELINES STILL UNDER CONSTRUCTION.

__Module 5b: full FoV, high resolution (sub-arcsec) imaging__

PIPELINES STILL UNDER CONSTRUCTION.












