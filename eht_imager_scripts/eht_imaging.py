#!/bin/python
import numpy as np
import ehtim as eh
import os

from ehtim.const_def import *
from ehtim.observing.obs_helpers import *
from ehtim.imaging.imager_utils import *

import argparse

def main(vis, closure_tels='DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609;PL610;PL611;PL612;IE613', npix=128, fov_arcsec=6., zbl=0., prior_fwhm_arcsec=1., doplots=False, imfile='myim', remove_tels='', niter=300, cfloor=0.001, conv_criteria=0.0001, use_bs=False ):

    ## for the generic pipeline, force the data types
    npix = int(npix)
    fov_arcsec  = float(fov_arcsec)
    zbl = float(zbl)
    prior_fwhm_arcsec = float(prior_fwhm_arcsec)
    niter = int(niter)
    cfloor = float(cfloor)
    conv_criteria = float(conv_criteria)

    vis1 = vis.rstrip('/') + '.ms' ## so the next line will work if it doesn't end in ms
    fitsout = vis1.replace('.MS','.fits').replace('.ms','.fits')

    ## remove any telescopes from the default list
    r_tels = remove_tels.split(';')
    for r_tel in r_tels:
	closure_tels = closure_tels.replace(r_tel,'')
    tmp_tel = closure_tels.split(';')
    tel = [ t for t in tmp_tel if t != '' ]

    ## starting from a measurement set, get a smaller ms with just the list of antennas
    ## use taql to get a list of telescopes
    ss = 'taql \'select NAME from %s/ANTENNA\' > closure_txt'%(vis)
    os.system(ss)
    os.system('grep -v select closure_txt > closure_which')
    idxtel = np.loadtxt('closure_which',dtype='S')
    atel = np.unique(np.ravel(tel))

    notfound = []
    aidx = np.array([],dtype='int')
    for a in atel:
        found_this = False
        for i in range(len(idxtel)):
            if a==idxtel[i][:len(a)]:
                aidx = np.append(aidx,i)
                found_this = True
        if not found_this:
            notfound.append (a)
    if len(notfound):
        print 'The following telescopes were not found:',notfound

    os.system( 'rm closure_*' )

    aidx_s = np.sort(aidx)

    if os.path.exists ('cl_temp.ms'):
        os.system('rm -fr cl_temp.ms')
    ss = 'taql \'select from %s where ' % vis
    for i in range (len(aidx_s)):
        for j in range (i+1, len(aidx_s)):
            ss += ('ANTENNA1==%d and ANTENNA2==%d' %(aidx_s[i],aidx_s[j]))
            if i==len(aidx_s)-2 and j==len(aidx_s)-1:
                ss += (' giving cl_temp.ms as plain\'')
            else:
                ss += (' or ')
    print 'Selecting smaller MS cl_temp.ms, this will take about 4s/Gb:'
    os.system (ss)

    ## check if weights are appropriate
    print 'Checking if weights are sensible values, updating them if not'
    ss = "taql 'select WEIGHT_SPECTRUM from cl_temp.ms limit 1' > weight_check.txt"
    os.system(ss)
    with open( 'weight_check.txt', 'r' ) as f:
        lines = f.readlines()
    f.close()
    first_weight = np.float(lines[3].lstrip('[').split(',')[0])
    if first_weight < 1e-9:
        ss = "taql 'update cl_temp.ms set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM*1e11'"
        os.system(ss)

    ## find the zero baseline amplitudes if it isn't set
    if zbl == 0:
	print 'Calculating the zero baseline flux (mean of amplitudes on DE601 -- DE605'
	## use baseline DE601 -- DE605 (shortest international to interational baseline)
	de601_idx = aidx[[ i for i, val in enumerate(tel) if 'DE601' in val ]][0]
	de605_idx = aidx[[ i for i, val in enumerate(tel) if 'DE605' in val ]][0]
	ss = "taql 'select means(gaggr(abs(DATA)),0,1) from cl_temp.ms' where ANTENNA1==%s and ANTENNA2=%s > amp_check.txt"%(de601_idx, de605_idx)
	os.system(ss)
	with open( 'amp_check.txt', 'r' ) as f:
	    lines = f.readlines()
        f.close()
        amps = np.asarray(lines[2].lstrip('[').rstrip(']\n').split(', '),dtype=float)
        zbl = np.max(amps)

    os.system('rm *_check.txt')

    ## convert vis to uv-fits
    ss = 'ms2uvfits in=cl_temp.ms out=%s writesyscal=F'%(fitsout)
    os.system(ss)

    ## observe the uv-fits
    print 'Reading in the observation to the EHT imager.'
    obs = eh.obsdata.load_uvfits(fitsout)

    if doplots:
        ## make plots
        print 'Plotting u-v coverage and amplitude vs. uv-distance'
        obs.plotall('u','v', conj=True, show=False, export_pdf=fitsout.replace('fits','u-v.pdf') ) ## u-v coverage
        obs.plotall('uvdist','amp', show=False, export_pdf=fitsout.replace('fits','uvdist-amp.pdf') ) ## etc

    ## the eht imager's default units is microarcseconds
    fov = fov_arcsec * 1e6 * RADPERUAS

    ## the dirty beam
    print 'Calculating the dirty beam.'
    dbeam = obs.dirtybeam(npix,fov)
    dbeam.save_fits(fitsout.replace('fits','dirty_beam.fits'))

    ## the clean beam
    print 'Calculating the clean beam.'
    cbeam = obs.cleanbeam(npix,fov)
    cbeam.save_fits(fitsout.replace('fits','clean_beam.fits'))

    # dirty image
    print 'Calculating the dirty image.'
    dimage = obs.dirtyimage(npix,fov)
    dimage.save_fits(fitsout.replace('fits','dirty_image.fits'))

    ## beam parameters (in radians)
    beamparams = obs.fit_beam()
    print 'The beam is ' + str(beamparams[0]*206265.) + ' by ' + str(beamparams[1]*206265.) + ' arcsec.'

    ## maximum resolution (in radians)
    res = obs.res()
    print 'Maximum resolution is ' + str(res*206265.)

    ## set up the gaussian prior
    prior_fwhm = prior_fwhm_arcsec * 1e6 * RADPERUAS
    gaussparams = ( prior_fwhm, prior_fwhm, 0.0 )
    emptyprior = eh.image.make_square( obs, npix, fov )
    gaussprior = emptyprior.add_gauss( zbl, gaussparams )
    gaussprior.display()

    if use_bs:
	print 'Using bispectrum mode rather than cphase and camp separately.'

	## try using bispectrum
        bs_out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='bs', s1='gs', alpha_d1=50, clipfloor=0.001, maxit=300, stop=0.0001 )
        bs_outblur = bs_out.blur_gauss(beamparams, 0.5)
        bs_out1 = eh.imager_func( obs, bs_outblur, bs_outblur, zbl, d1='bs', s1='gs', alpha_d1=50, clipfloor=0.001, maxit=300, stop=0.0001 )
        bs_outblur1 = bs_out1.blur_gauss((res,res,0),0.5)

        bs_finalout = bs_out1.blur_gauss(beamparams,0.5)

        ## save to fits
        bs_out1.save_fits('./bs_' + imfile + 'im.fits')
        bs_finalout.save_fits('./bs_' + imfile + 'im_blur.fits')

    else:
        print 'Using closure phase and closure amplitude separately.'
        ## using d1 = closure phase and d2 = closure amplitude
        out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50, alpha_d2=50, clipfloor=cfloor, maxit=niter, stop=conv_criteria )
        outblur = out.blur_gauss(beamparams, 0.5)
        out1 = eh.imager_func( obs, outblur, outblur, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50, alpha_d2=50, clipfloor=cfloor, maxit=niter, stop=conv_criteria )
        out1blur = out1.blur_gauss((res,res,0),0.5)

        finalout = out1.blur_gauss(beamparams,0.5)

        ## save to fits
	out.save_fits('./' + imfile + '_0_im.fits')
	outblur.save_fits('./' + imfile + '_0_im_blur.fits')

        out1.save_fits('./' + imfile + 'im.fits')
        finalout.save_fits('./' + imfile + 'im_blur.fits')

    print 'done.'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Python wrapper for running the EHT closure phase/closure amplitude imager on LOFAR data')
    parser.add_argument('vis', type=str, help='Measurement set name')
    parser.add_argument('--ctels', type=str, help='list of closure telescopes, separated by ; eg. "DE601;DE602" [default all international stations]', default='DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609;PL610;PL611;PL612;IE613')
    parser.add_argument('--npix', type=int, help='Number of pixels for output images', default=128)
    parser.add_argument('--fov', type=float, help='FoV in arcseconds', default=6. )
    parser.add_argument('--zbl', type=float, help='Zero baseline flux, Jy [default=calculate from data]', default=0. )
    parser.add_argument('--prior_fwhm', type=float, help='FWHM of Gaussian prior in arcsec', default=1. )
    parser.add_argument('--doplots', action='store_true')
    parser.add_argument('--imfile', type=str, help='filestem for output files', default='ehtim' )
    parser.add_argument('--rtels', type=str, help='list of telescopes to remove from ctels (same format)', default='')
    parser.add_argument('--maxiter', type=int, help='Maximum number of iterations', default=300)
    parser.add_argument('--clipfloor', type=float, help='Clip floor for model [Jy/pixel]', default=0.001)
    parser.add_argument('--convg', type=float, help='Convergence criteria', default=0.0001)
    parser.add_argument('--use_bs', action='store_true', help='set to use bs rather than cphase and camp separately')

    args = parser.parse_args()
    main( args.vis, closure_tels = args.ctels, npix = args.npix, fov_arcsec = args.fov, zbl = args.zbl, prior_fwhm_arcsec = args.prior_fwhm, doplots = args.doplots, imfile = args.imfile, remove_tels = args.rtels, niter=args.maxiter, cfloor=args.clipfloor, conv_criteria=args.convg, use_bs=args.use_bs )
