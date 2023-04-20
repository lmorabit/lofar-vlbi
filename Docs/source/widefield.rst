==================
Wide-field Imaging
==================

This is the step by step documentation for producing a wide-field image. 

#. Step 1.     The first step to producing a wide field image is to run LINC processing up to the end of LINC Target.
	
	*      The LINC initial calibration pipeline can be found here: https://git.astron.nl/RD/LINC 
	  
	*      For information on how to run the LINC pipeline please refer to the documentation: https://linc.readthedocs.io 

#. Step 2.     Now run the DDF pipeline which is needed for the later step of subtracting sources from outside the field. 

	*      This pipeline can be found here: https://github.com/mhardcastle/ddf-pipeline 

	*      You will need the appropriate Singularity image to run the pipeline inside of. To find the latest Singularity image go to https://lofar-webdav.grid.sara.nl/software/shub_mirror/tikk3r/lofar-grid-hpccloud/ and sort by last modified. Versions 4.0.0 + are compatible with cwl 

	*      Command to run pipeline: **pipeline.py  your_config_file.cfg** 

#. Step 3.     It’s time to start the Delay Calibration! Only steps up to applyddf DIS2 solutions and before running phaseup should be run. 

	*      The Delay Calibration parset can be found at https://github.com/lmorabit/lofar-vlbi 

	*      DIS2 solutions from running the DDF pipeline need to be applied so some alterations to the parset need to be made. Apply_ddf needs to be set along with the ddf_soldir which points to the solution ms files from the DDF run. Make sure that the DIS2 apply is on the 2MHz chunk ms files. It is also important to change the column from *DATA* to *DATA_DI_CORRECTED*.

	*      The pipeline steps also need to be updated so that phaseup is not run yet

#. Step 4.     In preparation for step 5, we need to make a box of a DS9 region we want to keep.  

	*      The python script make_box.py is available here: https://github.com/tikk3r/lofar-highres-widefield/blob/restructure/utils/make_box.py. This script should be run inside a Singularity image. 

	*      The script may need some slight adjustments. In case of errors try removing **ms=mslist** and hard code the path to a *msdpppconcat* file which was created in the Delay Calibration stage, so **ms=’path to target msdpppconcat files’**

	*      The recommended box size is a 2.5 degree x 2.5 degree region 

	*      To run the python script the command is: **Python make_box.py $path_to_target_msdpppconcat_files$ 2.5**

#. Step 5.     The next stage is the DDF source subtraction step. 

	*      The code can be found here: https://github.com/rvweeren/lofar_facet_selfcal/blob/main/sub-sources-outside-region.py and should be run within the singularity image 

	*      Several files/directories are needed from the DDF pipeline run: *SOLSDIR*, *DDS3*.npz*, *image_dirin_SSD_m.npy.ClusterCat.npy*, *image_full_ampphase_di_m.NS.DicoModel*, *image_full_ampphase_di_m.NS.mask01.fits*,  *big-mslist.txt* and *mslist.txt*    

	*      The files that the script works on are the *pre-cal.ms* directory from *SOLSDIR* and the *msdpppconcat* measurement sets from the Delay Calibrator run before phaseup is ran (step 3). *SOLSDIR* and the *msdpppconcat* directories should be within the DDF working directory. 

	*      It is recommended to make copies of the *msdpppconcat* files as the **sub-sources-outside-region.py** alters these measurement sets 

	*      The names of the *msdpppconcat* directories must match EXACTLY the names of the files in *SOLSDIR*. So, the names of *msdppconcat* have to be altered accordingly. The file names should then be inputted into the *mslist.txt* file. Each subtraction of a single 2MHz chunk can take up to 20 hours so it is recommended to use a cluster where possible to run several chunks of different nodes, so you may want only one subband in *mslist.txt*. 

	*      Command to run in singularity: **python sub-sources-outside-region.py --boxfile boxfile.reg --column DATA_DI_CORRECTED --freqavg 1 --timeavg 1 --ncpu 24 --prefixname sub6asec --noconcat --keeplongbaseline --nophaseshift --chunkhours 1 --onlyuseweightspectrum --mslist mslist.txt** 

	*      If you encounter memory issues lower the value of chunkhours and rerun

#. Step 6.     Verify that the subtraction was successful by imaging three blocks across the band or all blocks. 

	*      These should match the residual of ddf-pipeline 

#. Step 7.     Now complete the rest of delay calibration from phase up step by removing **#** from steps to run in second figure 

#. Step 8.     Check delay calibration solutions and apply within the pipeline if reasonable. 

#. Step 9.     A 1” image can now be made. This step is important to check the quality of DI solutions. This is done using one wsclean command , below is an example:

	*      **wsclean -update-model-required -minuv-l 80.0 -size 22500 22500 -weighting-rank-filter 3 -reorder -weight briggs -1.5 -parallel-reordering 6 -mgain 0.65 -data-column DATA -auto-mask 3 -auto-threshold 1.0 -pol i -name 1.2asec_I -scal0.4arcsec -taper-gaussian 1.2asec -niter 50000 -log-time -multiscale-scale-bias 0.6 -parallel-deconvolution 2600 -multiscale -multiscale-max-scales 9 -nmiter 1 -mem 25 -channels-out 1 -j 32 -use-idg -grid-wit-beam -use-differential-lofar-beam <input ms> **

#. Step 10.    Split out DDE calibrators  

#. Step 11.    Using **facetselfcal** we now do self-calibration on DDE calibrators  

#. Step 12.    Remake 1” image with DDE calibration 

#. Step 13.     The final step is to make the 0.3” image. This step very much depends on compute power available 
