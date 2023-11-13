

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

What it does: Queries 













