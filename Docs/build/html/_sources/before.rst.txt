.. index:: Before you begin

****************
Before you begin
****************

=================
Selecting a field 
=================

Before beginning any data processing, it is important to keep several things in mind. First, the pipeline assumes that you are giving it something sensible, and it does not (yet) make any quality checks. That means if you give it data which it cannot handle, you will not get quality results. Before you select a field, please follow these guidelines:

* **Ionosphere:** should not be particularly bad. Check inspection plots from the observatory.  If it's been processed successfully by the Surveys KSP infrastructure, it should be okay.
* **Field of View:** The Surveys standard averaging is 16 channels/SB and 1 sec, which limits the Field of View (FoV) to ~1.25 deg from the centre. Your target and LBCS calibrators should be inside this.
* **Distance between calibrator and target:** An acceptable value will vary based on ionosphere and the distribution of sources in the field. As a rule of thumb, 1 degree separation *should* be okay, but it is safer to look for a calibrator within 0.5 degrees. The pipeline will pick the best in-field calibrator if there are multiple choices, but you can always manually phase reference a calibrator closer to your target from a more distant calibrator.
* **Flux Density:** In principle the limiting rms noise will be ~100 uJy/bm for an 8 hour observation. In practice, the quality of your image / amount you have to self-calibrate will depend on all of the factors above, including the distance from the phase centre of the observation. 
* **Source structure:** Keep in mind that not all point sources at 6" resolution are point sources at 0.3" resolution. If you expect your 1 mJy/bm source to break equally into 4 components, each component will only be 250 uJy/bm and therefore only a 2.5 sigma detection for 100 uJy/bm rms. 
* **Flux calibrator:** 3C295 and 3C196 are commonly used as standard flux calibrators, but high-resolution models of them are not included in the pipeline. They are therefore currently not advisable to use as flux calibrators for LOFAR-VLBI.

As an example, a good field will have: the target within ~ 1 degree from the pointing centre, an LBCS calibrator < 0.5 degrees from your target, no bright (> 1 Jy) sources within ~ 1 degree, reasonable to good ionospheric conditions, with 16 channels/SB and 1 second time resolution. If your data is averaged more than this, you may run into problems. 

=============================
Helpful information and links
=============================

The `LOFAR-VLBI`_ pipeline, while self-contained, works within a larger framework of LOFAR data processing. 
As the LOFAR software installation and maintenance can be quite complex, we also **strongly recommend** the use of `Singularity`_ for your data processing. Before you begin, please take some time to familiarise yourself with the following:

* **Singularity** is the way we *containerize* the LOFAR software. That means all the LOFAR software you need is included in one Singularity 'image' which you invoke to run the pipeline. You do not need to build your own image, but it is useful to understand how Singularity works: see `Singularity`_ .

* **Prefactor** is used for direction-independent calibration of the Dutch stations, prior to running the LOFAR-VLBI pipeline. It is useful to familiarise yourself with `prefactor`_ before using it. A `prefactor tutorial`_ and the `prefactor documentation`_ are helpful places to start. 

The `LOFAR-VLBI`_ pipeline makes use of two catalogue servers. It is a good idea to check whether or not there is coverage of your field beforehand:

* **LBCS:** The Long Baseline Calibrator Survey (LBCS) is now complete, and provides information on LOFAR-VLBI in-field calibrators across the sky. You should check that there is coverage for your field before you process your data. The catalogue can be searched here: `LBCS catalogue server`_

* **LoTSS:** Useful information is collected from the LOFAR Two-metre Sky Survey (LoTSS) on the sources in your field. If LoTSS does not cover your field, you will have to manually produce a catalogue at 6 arcsecond resolution using the Dutch array before you run the LOFAR-VLBI pipeline. You can check the LoTSS catalogue here: `LoTSS catalogue server`_

   
.. _help:

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
