#!/usr/bin/env python
import os
import argparse
import pyrap.tables as pt
import numpy as np

def writePhaseOnlyParset( parsetname="phaseonly.parset", incol="DATA", outcol="CORRECTED_DATA", solint=10, outms="." ):
    # open write-only parset; this will overwrite previous files of the same name
    with open(parsetname,"w") as f:
        f.write("msin.datacolumn=%s\n"%incol)
        f.write("msout=%s\n"%outms)
	f.write("msout.datacolumn=%s\n"%outcol)
        f.write("steps=[gaincal]\n")
        f.write("gaincal.usemodelcolumn=true\n")
        f.write("gaincal.caltype=phaseonly\n")
        f.write("gaincal.solint=%i\n"%(solint))
        f.write("gaincal.maxiter=200\n")
        f.write("gaincal.tolerance=0.0001\n")
	f.write("gaincal.applysolution=True\n")
    f.close()

def writePhaseAmpParset( parsetname="phaseamp.parset", incol="DATA", outcol="CORRECTED_DATA", solint=10, outms="." ):
    # open write-only parset; this will overwrite previous files of the same name
    with open(parsetname,"w") as f:
        f.write("msin.datacolumn=%s\n"%incol)
        f.write("msout=%s\n"%outms)
	f.write("msout.datacolumn=%s\n"%outcol)
        f.write("steps=[gaincal]\n")
        f.write("gaincal.usemodelcolumn=true\n")
        f.write("gaincal.caltype=diagonal\n")
        f.write("gaincal.solint=%i\n"%(solint))
        f.write("gaincal.maxiter=200\n")
        f.write("gaincal.tolerance=0.0001\n")
	f.write("gaincal.applysolution=True\n")
    f.close()

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'msname', type=str, help='Measurement set name' )
    #parser.add_argument( '-m', '--msname', dest='msname', type=str, help='Measurement set name', required=True )
    parser.add_argument( '-l', dest='large_im', action='store_true', default=False, help='Turn on to make a larger image of the field.' )
    parser.add_argument( '-f', '--freq_range', dest='freq_range', default=10, help='Frequency range in MHz before breaking into channels to image (default 10)' )
    parser.add_argument( '-p', dest='do_phaseonly', action='store_true', default=False, help='Run phase only self-calibration' )
    parser.add_argument( '-a', dest='do_ampphase', action='store_true', default=False, help='Run amplitude and phase only self-calibration' )
    parser.add_argument( '-d', dest='data_col', type=str, default='DATA', help='Which data column to start with the imaging' )

    args = parser.parse_args()
    msname = args.msname
    large_im = args.large_im
    freq_range = args.freq_range
    do_phaseonly = args.do_phaseonly
    do_ampphase = args.do_ampphase
    data_col = args.data_col

    ## get the target name
    tgtname = msname.split('.ms')[0].split('_')[0]

    ## find the frequency range to determine if MFS needs to be used
    t = pt.table( msname.rstrip('/') + '/SPECTRAL_WINDOW' )
    bandwidth = t.getcol('TOTAL_BANDWIDTH')[0] / 1e6 ## convert to MHz
    ## check if the bandwidth is smaller or larger than the allowed freq_range
    if bandwidth > freq_range:
        ## split into channels
        nchan = int( np.ceil( bandwidth / freq_range ) )

    ## INITIAL DIRTY IMAGES AT THREE DIFFERENT SCALES
    print "Inverting visibilities at 3 scales to check for signal."
    ## wsclean can take a list of measurement sets -- they don't need to be concatenated prior to running wsclean
    ## 800L = 1.6 km @ 150 MHz - theta = 260 arcsec -- this is not worthwhile imaging, since we combine core stations.
    ## 8000L = 16 km @ 150 MHz - theta = 26 arcsec
    ## 80000 = 159 km @ 150 MHz - theta = 2.6 arcsec
    ## 800000 = 1590 km @ 150 MHz - theta = 0.26 arcsec

    pix_scales = [ '5asec', '0.5asec', '0.05asec' ]
    max_lambda = [ 8000, 80000, 800000 ]
    sizes = [ 128, 1024, 4096 ]

    for ii in range(0,3):
        ## set the parameters
	    scale = pix_scales[ii]
	    maxuv = max_lambda[ii]
	    imsize = sizes[ii]*2
	    imtrim = sizes[ii]
        ## run wsclean
	    ss = "wsclean -datacolumn DATA -make-psf -scale %s -size %i %i -trim %i %i -maxuv-l %i -reorder -name %s-%sL-NOCLN %s"%(scale,imsize,imsize,imtrim,imtrim,maxuv,tgtname,str(maxuv),msname)
	    os.system( ss )

    ## FIRST CLEAN IMAGE
    print "Starting first clean."
    ## set the parameters
    if large_im:
        imtrim = 7776
    else:
        imtrim = 4096
    ## get the right size
    imsize = imtrim*2

    ## check if MFS needs to be turned on
    if bandwidth > freq_range:
        ss = "wsclean -datacolumn %s -scale 0.05asec -size %i %i -trim %i %i -joinchannels -channelsout %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s-CLN %s"%(data_col,imsize,imsize,imtrim,imtrim,nchan,tgtname,msname)
    else:
        ss = "wsclean -datacolumn %s -scale 0.05asec -size %i %i -trim %i %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s-CLN %s"%(data_col,imsize,imsize,imtrim,imtrim,tgtname,msname)

    ## run imaging
    os.system( ss )
 
    if do_phaseonly:
    ## SET UP FOR PHASE-ONLY SELF-CALIBRATION
        ## remove previous log
        os.system( 'rm ndppp_%s.log'%(tgtname) )

        ## ranges to use
        bl_min = [ 0, 0 ]
        bl_max = [ 1e30, 1e30 ]
        solints = [ 4, 2 ]

        nloops = len( bl_max )

        for ii in range(nloops):
    
            print "Performing phase-only calibration " + str( ii+1 )
            ## probably want to change the bl-length and solint for more iterations
            writePhaseOnlyParset( solint=solints[ii] )
            ss = "NDPPP phaseonly.parset msin=%s gaincal.parmdb=%s/phaseonly-%s gaincal.blrange=[%i,%i] >> ndppp_%s.log 2>&1"%(msname,msname,str(ii+1),bl_min[ii],bl_max[ii],tgtname)
            os.system( ss )
            imagename = "%s-CLN-PO%s"%(tgtname,str(ii+1))
            if bandwidth > freq_range:
                ss = "wsclean -datacolumn CORRECTED_DATA -scale 0.05asec -size %i %i -trim %i %i -joinchannels -channelsout %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s-CLN %s"%(imsize,imsize,imtrim,imtrim,nchan,tgtname,msname)
            else:
                ss = "wsclean -datacolumn CORRECTED_DATA -scale 0.05asec -size %i %i -trim %i %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s %s"%(imsize,imsize,imtrim,imtrim,imagename,msname)
            os.system( ss )

    if do_ampphase:

        ## PHASE AND AMP
        writePhaseAmpParset()
        ss = "NDPPP phaseamp.parset msin=%s gaincal.parmdb=%s/phaseamp gaincal.blrange=[%i,%i] >> ndppp_%s.log 2>&1"%(msname,msname,bl_min[-1],bl_max[-1],tgtname)
        os.system( ss )
        imagename = "%s-CLN-PA-final"%(tgtname)
        if bandwidth > freq_range:
            ss = "wsclean -datacolumn CORRECTED_DATA -scale 0.05asec -size %i %i -trim %i %i -joinchannels -channelsout %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s-CLN %s"%(imsize,imsize,imtrim,imtrim,nchan,tgtname,msname)
        else:
            ss = "wsclean -datacolumn CORRECTED_DATA -scale 0.05asec -size %i %i -trim %i %i -rms-background -maxuv-l 800000 -gain 0.1 -mgain 0.85 -niter 100000 -weighting-rank-filter 3 -weighting-rank-filter-size 256 -auto-threshold 5 -multiscale -auto-mask 7 -save-component-list -reorder -name %s %s"%(imsize,imsize,imtrim,imtrim,imagename,msname)
        os.system( ss )

if __name__=="__main__":
    main()

