#!/bin/python
import numpy as np
import ehtim as eh
import os
import sys

from ehtim.const_def import *
from ehtim.observing.obs_helpers import *
from ehtim.imaging.imager_utils import *


def main(vis, closure_tels='DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609', npix=128, fov_arcsec=6., zbl=350., prior_fwhm_arcsec=1., doplots=False, imfile='3C48' )

    fitsout = vis.replace('MS','fits').replace('ms','fits')

    ## starting from a measurement set, get a smaller ms with just the list of antennas
    closure_tels = 'DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609'
    tel = closure_tels.split(';')

    ## use taql to get a list of telescopes
    ss = 'taql \'select NAME from %s/ANTENNA\' >closure_txt'%(vis)
    os.system(ss)
    os.system('grep -v select closure_txt >closure_which')
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
        ss = "taql 'select WEIGHT_SPECTRUM from cl_temp.ms limit 1' > weight_check.txt"
        os.system(ss)
        with open( 'weight_check.txt', 'r' ) as f:
            lines = f.readlines()
        f.close()
        first_weight = np.float(lines[3].lstrip('[').split(',')[0])
    if first_weight < 1e-9:
        ss = "taql 'update cl_temp.ms set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM*1e11'"
        os.system(ss)

    ## convert vis to uv-fits
    ss = 'ms2uvfits in=cl_temp.ms out=%s writesyscal=F'%(fitsout)
    os.system(ss)

    ## observe the uv-fits
    obs = eh.obsdata.load_uvfits(fitsout)

    if doplots:
        ## make plots
        obs.plotall('u','v', conj=True, show=False, export_pdf=fitsout.replace('fits','u-v.pdf') ) ## u-v coverage
        obs.plotall('uvdist','amp', show=False, export_pdf=fitsout.replace('fits','uvdist-amp.pdf') ) ## etc

    ## the eht imager's default units is microarcseconds
    fov = fov_arcsec * 1e6 * RADPERUAS

    ## the dirty beam
    dbeam = obs.dirtybeam(npix,fov)
    dbeam.save_fits(fitsout.replace('fits','dirty_beam.fits'))

    ## the clean beam
    cbeam = obs.cleanbeam(npix,fov)
    cbeam.save_fits(fitsout.replace('fits','clean_beam.fits'))

    # dirty image
    dimage = obs.dirtyimage(npix,fov)
    dimage.save_fits(fitsout.replace('fits','dirty_image.fits'))

    ## beam parameters (in radians)
    beamparams = obs.fit_beam()

    ## maximum resolution (in radians)
    res = obs.res()

    ## set up the gaussian prior
    prior_fwhm = prior_fwhm_arcsec * 1e6 * RADPERUAS
    gaussparams = ( prior_fwhm, prior_fwhm, 0.0 )
    emptyprior = eh.image.make_square( obs, npix, fov )
    gaussprior = emptyprior.add_gauss( zbl, gaussparams )
    gaussprior.display()

    ## using d1 = closure phase and d2 = closure amplitude
    out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50, alpha_d2=50, clipfloor=0.001, maxit=300, stop=0.0001 )
    outblur = out.blur_gauss(beamparams, 0.5)
    out1 = eh.imager_func( obs, outblur, outblur, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50, alpha_d2=50, clipfloor=0.001, maxit=300, stop=0.0001 )
    out1blur = out1.blur_gauss((res,res,0),0.5)

    finalout = out1.blur_gauss(beamparams,0.5)

    ## save to fits
    out1.save_fits('./' + imfile + 'im.fits')
    finalout.save_fits('./' + imfile + 'im_blur.fits')

    ## try using bispectrum
    bs_out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='bs', s1='gs', alpha_d1=50, clipfloor=0.001, maxit=300, stop=0.0001 )
    bs_outblur = bs_out.blur_gauss(beamparams, 0.5)
    bs_out1 = eh.imager_func( obs, bs_outblur, bs_outblur, zbl, d1='bs', s1='gs', alpha_d1=50, clipfloor=0.001, maxit=300, stop=0.0001 )
    bs_outblur1 = bs_out1.blur_gauss((res,res,0),0.5)

    bs_finalout = bs_out1.blur_gauss(beamparams,0.5)

    ## save to fits
    bs_out1.save_fits('./bs_' + imfile + 'im.fits')
    bs_finalout.save_fits('./bs_' + imfile + 'im_blur.fits')


if __name__ == "__main__":

    myvis = sys.argv[1]
    main(vis)

