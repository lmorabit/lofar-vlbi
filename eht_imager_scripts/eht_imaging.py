#!/bin/python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import ehtim as eh
import os
import math 
import casacore.tables as casatb
import pyfits 
import argparse 
import aplpy
from matplotlib import pyplot as plt
from ehtim.const_def import *
from ehtim.observing.obs_helpers import *
from ehtim.imaging.imager_utils import *
from astropy.coordinates import SkyCoord
from astropy.io import ascii

# Plot a fitsfile giving a png out. Optionally specify the white/black levels. Return the 0.25% and
# 99.75% percentile; these are the aplpy defaults for plotting, so this allows a subsequent call to
# use the same colour scale

def gcmake (infits, outpng, vmin=None, vmax=None):
    a=pyfits.getdata(infits)
    gc = aplpy.FITSFigure(infits)
    gc.ticks.show()
    gc.show_colorscale(cmap=matplotlib.cm.Oranges,vmin=vmin,vmax=vmax)
    gc.set_title(infits)
    gc.axis_labels.set_font(size=18)
    gc.tick_labels.set_font(size=18)
    gc.add_label (0.5,0.05,'Peak %4f Jy'%a.max(),relative=True,size=14)
    gc.savefig(outpng)
    aa = np.sort(np.ravel(a))
    return aa[int(len(aa)*0.25/100)],aa[int(len(aa)*99.75/100)] # default aplpy scaling

# The EHT imager writes fitsfiles without CRVAL keywords: insert them

def insert_radec(fits,ra,dec):
    a,h = pyfits.getdata(fits),pyfits.getheader(fits)
    os.system('rm '+fits)
    h['CRVAL1'], h['CRVAL2'] = ra,dec
    pyfits.writeto(fits,a,h)

# Conversion rad<-> arcsec (makes subsequent code a lot tidier)

def as2rad (z):
    return z/(3600.*180/np.pi)

def rad2as (z):
    return z*3600.*180/np.pi

# FITS files delivered from FIRST have FREQ and STOKES axes which plotting software (in particular
# aplpy) is unable to cope with, so remove them.

def wcsreduce(fits):
    a,h=pyfits.getdata(fits),pyfits.getheader(fits)
    h.remove('CTYPE3');h.remove('CRVAL3');h.remove('CRPIX3');h.remove('CDELT3');h.remove('CROTA3')
    h.remove('CTYPE4');h.remove('CRVAL4');h.remove('CRPIX4');h.remove('CDELT4');h.remove('CROTA4')
    os.system('rm '+fits)   # do not use 'clobber' - non-backwards compatible change in pyfits
    pyfits.writeto(fits,a,h)
    

def first_download (ra,dec,outfile='first_out.fits',gif=0,fits=1,imsize=2.0,imserver='third.ucllnl.org'):
    os.system('rm '+outfile)
    coord = SkyCoord(ra,dec,unit='deg')
    scoord = str(coord.to_string(style='hmsdms'))
    s = scoord.replace('h',' ').replace('d',' ').replace('m',' ').replace('s',' ').split()
    command = ('wget -O %s "https://%s/cgi-bin/firstimage?RA=%s %s %s %s %s %s&Dec=&Equinox=J2000&ImageSize=%.1f&MaxInt=10&GIF=%d&FITS=%d&Download=1"'%(outfile,imserver,s[0],s[1],s[2],s[3],s[4],s[5],imsize,gif,fits))
    os.system(command)
    wcsreduce(outfile) # aplpy can't cope with FREQ,STOKES axes
    return pyfits.getdata(outfile)

# Altered to use astropy.SkyCoord.separation

def sepn(r1,d1,r2,d2):
    return SkyCoord(r1,d1,unit='rad').separation(SkyCoord(r2,d2,unit='rad')).rad


def main(vis, closure_tels='DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609;PL610;PL611;PL612;IE613', npix=128, fov_arcsec=6., maxarcmin=2.0, zbl=0., prior_fwhm_arcsec=1., doplots=False, imfile='myim', remove_tels='', niter=300, cfloor=0.001, conv_criteria=0.0001, use_bs=False, scratch_model=False, do_diagnostic_plots=False, lotss_file=''):

    use_first = not(scratch_model)

    ## for the generic pipeline, force the data types

    npix = int(npix)
    fov_arcsec  = float(fov_arcsec)
    maxarcmin = float(maxarcmin)
    zbl = float(zbl)
    prior_fwhm_arcsec = float(prior_fwhm_arcsec)
    niter = int(niter)
    cfloor = float(cfloor)
    conv_criteria = float(conv_criteria)

    ## set output files 

    vis1 = vis.rstrip('/') + '.ms' ## so the next line will work if it doesn't end in ms
    fitsout = vis1.replace('.MS','.fits').replace('.ms','.fits')
    tmp_name = 'tmp_clphase_' + vis1
    src_name = vis1.split('_')[0]
    
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
    os.system('cat closure_which')
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
    if os.path.exists (tmp_name):
        os.system('rm -fr %s'%tmp_name)
    ss = 'taql \'select from %s where ' % vis
    for i in range (len(aidx_s)):
        for j in range (i+1, len(aidx_s)):
            ss += ('ANTENNA1==%d and ANTENNA2==%d' %(aidx_s[i],aidx_s[j]))
            if i==len(aidx_s)-2 and j==len(aidx_s)-1:
                ss += (' giving %s as plain\''%tmp_name)
            else:
                ss += (' or ')
    print 'Selecting smaller MS %s this will take about 4s/Gb:'%tmp_name
    os.system (ss)

    ## get the coordinates for the phase center
    ss = "taql 'select PHASE_DIR from %s/FIELD' > phase_dir.txt"%tmp_name
    os.system( ss )
    with open( 'phase_dir.txt', 'r' ) as f:
	lines = f.readlines()
    f.close()

    phase_dir = lines[-1].lstrip('[').rstrip(']\n')
    ra_sexg = phase_dir.split(',')[0]
    dec_sexg = phase_dir.split(',')[1]
    c = SkyCoord( ra_sexg, dec_sexg, frame='icrs' )
    ra = float(c.ra.degree)
    dec = float(c.dec.degree)

    ## check if weights are appropriate

    print 'Checking if weights are sensible values, updating them if not'
    ss = "taql 'select WEIGHT_SPECTRUM from %s limit 1' > weight_check.txt"%tmp_name
    print( ss )
    os.system(ss)
    with open( 'weight_check.txt', 'r' ) as f:
        lines = f.readlines()
    f.close()
    first_weight = np.float(lines[3].lstrip('[').split(',')[0])
    if first_weight < 1e-9:
        ss = "taql 'update %s set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM*1e11'"%tmp_name
        os.system(ss)

    if lotss_file != '':
	print 'LoTSS file is specified, using zero baseline flux from catalogue.'
	lotss_cat = ascii.read( lotss_file )
	src_index = np.where( lotss_cat['Source_id'] == src_name )
	flux_mJy = float( lotss_cat['Total_flux'][src_index[0]] )
	zbl = flux_mJy * 1e-3
	
    else:
	print 'LoTSS file is not specified.'
	if zbl == 0:
            ## find the zero baseline amplitudes if it isn't set
  	    print 'Calculating the zero baseline flux (mean of amplitudes on DE601 -- DE605)'
	    ## use baseline DE601 -- DE605, failing that the first two in list
            try: 
	        de601_idx = aidx[[ i for i, val in enumerate(tel) if 'DE601' in val ]][0]
	        de605_idx = aidx[[ i for i, val in enumerate(tel) if 'DE605' in val ]][0]
	    except:
                print '...Not present. Using %s %s instead'%(tel[0],tel[1])
                de601_idx, de605_idx = aidx[0],aidx[1]
	    ss = "taql 'select medians(gaggr(abs(DATA)),0,1) from %s' where ANTENNA1==%s and ANTENNA2=%s > amp_check.txt"%(tmp_name, de601_idx, de605_idx)
	    os.system(ss)
	    with open( 'amp_check.txt', 'r' ) as f:
	        lines = f.readlines()
            f.close()
            amps = np.asarray(lines[2].lstrip('[').rstrip(']\n').split(', '),dtype=float)
            zbl = np.median(amps)
    	    os.system('rm *_check.txt')

    print( 'Zero baseline flux:', zbl )

    ## convert vis to uv-fits

    if os.path.isfile(fitsout):
	os.system('rm '+fitsout)
    ss = 'ms2uvfits in=%s out=%s writesyscal=F'%(tmp_name, fitsout)
    os.system(ss)
    ## observe the uv-fits

    print 'Reading in the observation %s to the EHT imager.' %fitsout
    obs = eh.obsdata.load_uvfits(fitsout)

    ## flag the data?
    fov = as2rad(fov_arcsec)

    if do_diagnostic_plots:
        print 'Plotting u-v coverage and amplitude vs. uv-distance'
        obs.plotall('u','v', conj=True, show=False, export_pdf=fitsout.replace('fits','u-v.pdf') ) ## u-v coverage
        obs.plotall('uvdist','amp', show=False, export_pdf=fitsout.replace('fits','uvdist-amp.pdf') ) ## etc

        ## the dirty beam
	print 'Calculating the dirty beam.'
        dbeam = obs.dirtybeam(npix,fov)
        dbeam.save_fits(fitsout.replace('fits','dirty_beam.fits'))

        ## the clean beam
        print 'Calculating the clean beam.'
        cbeam = obs.cleanbeam(npix,fov)
        cbeam.save_fits(fitsout.replace('fits','clean_beam.fits'))

        ## dirty image
        print 'Calculating the dirty image.'
        dimage = obs.dirtyimage(npix,fov)
        dimage.save_fits(fitsout.replace('fits','dirty_image.fits'))

    ## beam parameters (in radians)
    beamparams = obs.fit_beam()
    print 'The beam is ' + str(rad2as(beamparams[0])) + ' by ' + str(rad2as(beamparams[1])) + ' arcsec.'

    ## maximum resolution (in radians)
    res = obs.res()
    print 'Maximum resolution is ' + str(res*206265.)
    print use_first
    ## set up the gaussian prior
    if use_first:
	print( 'Using FIRST.' )
        if not os.path.isfile('./first_2008.simple.npy'):
            os.system('wget http://www.jb.man.ac.uk/~njj/first_2008.simple.npy')
        first = np.load('first_2008.simple.npy')
        #firstdata = first_download(ra,dec,imsize=maxarcmin)
	# lkm - this is the else statement, it may not work 
	spatial_separations = sepn(np.deg2rad(ra),np.deg2rad(dec),np.deg2rad(first[:,0]),np.deg2rad(first[:,1]))
	corrfirst = np.where( spatial_separations <= as2rad(maxarcmin*60.0) )[0]
        if not len(corrfirst):           # no FIRST, fall back to default
            print 'No FIRST sources, using default'
            prior_fwhm = as2rad(prior_fwhm_arcsec)
            gaussparams = ( prior_fwhm, prior_fwhm, 0.0 )
            emptyprior = eh.image.make_square( obs, npix, fov )
            gaussprior = emptyprior.add_gauss( zbl, gaussparams )
            print 'Prior image: circular Gaussian %f arcsec' % prior_fwhm_arcsec
        else:
            gaussparams, gflux = np.array([]),np.array([])
            for i in corrfirst:
                this = first[int(i)]
                gflux = np.append(gflux,this[3])
                decdiff = np.deg2rad(this[1]-dec)
                radiff = np.deg2rad(this[0]-ra)*np.cos(np.deg2rad(dec))
                gnew = np.array([as2rad(max(7.0,this[4])), as2rad(max(5.0,this[5])),\
                                 np.deg2rad(this[6]),radiff,decdiff])
                try:
                    gaussparams = np.vstack((gaussparams,gnew))
                except:
                    gaussparams = np.copy([gnew])
                fov_arcsec = max(fov_arcsec,2.5*np.hypot(rad2as(gnew[3]),rad2as(gnew[4])))
            print 'Adjusted FOV to %f arcsec'%fov_arcsec
#            gaussparams[0,4]*=0.8; gaussparams[0,3]*=0.8       #sabotage
            fov = as2rad(fov_arcsec)
            gaussprior = eh.image.make_square( obs, npix, fov )
            zbl = np.sum(gflux) * 4.0/1000  # total flux in FIRST * 4
            for i in range(len(gaussparams)):
                g = gaussparams[i]
                print 'Prior: adding %.1fmJy Gaussian %.2f*%.2f asec, PA %.1f, at (%.3f,%.3f)asec' % \
                  (gflux[i],rad2as(g[0]),rad2as(g[1]),np.rad2deg(g[2]),rad2as(g[3]),rad2as(g[4]))
                gaussprior = gaussprior.add_gauss( gflux[i], g )
#            firstdata = first_download(ra,dec,imsize=fov_arcsec/60.0)
    else:
        print 'Not using FIRST'
#        obs=eh.obsdata.load_uvfits('../data/L1327+5504_vvsmall.fits') ##
        prior_fwhm = as2rad(prior_fwhm_arcsec)
        gaussparams = ( prior_fwhm, prior_fwhm, 0.0 )
        emptyprior = eh.image.make_square( obs, npix, fov )
        print npix,fov,zbl,gaussparams
        gaussprior = emptyprior.add_gauss( zbl, gaussparams )
        out=eh.imager_func(obs,gaussprior,gaussprior,5.0,d1='bs',clipfloor=cfloor)


    if use_bs:
	print 'Using bispectrum mode rather than cphase and camp separately.'
        out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='bs', \
                              clipfloor=cfloor, maxit=300)
        outblur = out.blur_gauss(beamparams, 0.5)
        out1 = eh.imager_func( obs, outblur, outblur, zbl, d1='bs', \
                            clipfloor=cfloor, maxit=300 )
#        out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='bs', s1='gs', \
#              alpha_d1=50, clipfloor=0.00, maxit=300, stop=0.0001 )
#       outblur = out.blur_gauss(beamparams, 0.5)
#        out1 = eh.imager_func( obs, outblur, outblur, zbl, d1='bs', s1='gs', \
#               alpha_d1=50, clipfloor=0.00, maxit=300, stop=0.0001 )
    else:
        print 'Using closure phase and closure amplitude separately.'
        ## using d1 = closure phase and d2 = closure amplitude
        out = eh.imager_func( obs, gaussprior, gaussprior, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50, alpha_d2=50, clipfloor=cfloor, maxit=niter, stop=conv_criteria )
        outblur = out.blur_gauss(beamparams, 0.5)
        out1 = eh.imager_func( obs, outblur, outblur, zbl, d1='cphase', d2='camp', s1='gs', s2='gs', alpha_d1=50,alpha_d2=50, clipfloor=cfloor, maxit=niter, stop=conv_criteria )

    outblur1 = out1.blur_gauss((res,res,0),0.5)
    finalout = out1.blur_gauss(beamparams,0.5)

# saving final fits file and some intermediate ones

    gaussprior.save_fits('./' + imfile + 'prior.fits')
    out.save_fits('./' + imfile + '_0_im.fits')
    outblur.save_fits('./' + imfile + '_0_im_blur.fits')
    out1.save_fits('./' + imfile + 'im.fits')
    finalout.save_fits('./' + imfile + 'im_blur.fits')
    insert_radec ('./' + imfile + 'prior.fits', ra, dec)
    insert_radec ('./' + imfile + 'im_blur.fits', ra, dec)

    if doplots and use_first:
    #    a1,a2 = gcmake ('first_out.fits','./%s_p1.png'%imfile)
        a1,a2 = gcmake ('./%sprior.fits'%imfile,'./%s_p2.png'%imfile)
        a1,a2 = gcmake ('./%sim_blur.fits'%imfile,'./%s_p3.png'%imfile,vmax=a2)
    #    os.system('montage -geometry 600x600 -tile 3x1 %s_p1.png %s_p2.png %s_p3.png %s_plots.png'%(imfile,imfile,imfile,imfile))
    if doplots and not use_first:
            plt.subplot(121);plt.imshow(pyfits.getdata('./'+imfile+'im_blur.fits'))
            plt.subplot(122);plt.imshow(pyfits.getdata('./'+imfile+'prior.fits'))
            plt.savefig('./'+imfile+'_plots.png',bbox_inches='tight')
    print 'done.'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Python wrapper for running the EHT closure phase/closure amplitude imager on LOFAR data')
    parser.add_argument('vis', type=str, help='Measurement set name')
    parser.add_argument('--ctels', type=str, help='list of closure telescopes, separated by ; eg. "DE601;DE602" [default all international stations]', default='DE601;DE602;DE603;DE604;DE605;FR606;SE607;UK608;DE609;PL610;PL611;PL612;IE613')
    parser.add_argument('--npix', type=int, help='Number of pixels for output images', default=128)
    parser.add_argument('--fov', type=float, help='FoV in arcseconds', default=6. )
    parser.add_argument('--maxarcmin', type=float, help='Maximum arcmin size of field', default=2.0 )
    parser.add_argument('--zbl', type=float, help='Zero baseline flux, Jy [default=calculate from data]', default=0. )
    parser.add_argument('--prior_fwhm', type=float, help='FWHM of Gaussian prior in arcsec', default=1. )
    parser.add_argument('--doplots', action='store_true')
    parser.add_argument('--diagplots', action='store_true')
    parser.add_argument('--imfile', type=str, help='filestem for output files', default='ehtim' )
    parser.add_argument('--rtels', type=str, help='list of telescopes to remove from ctels (same format)', default='')
    parser.add_argument('--maxiter', type=int, help='Maximum number of iterations', default=300)
    parser.add_argument('--clipfloor', type=float, help='Clip floor for model [Jy/pixel]', default=0.001)
    parser.add_argument('--convg', type=float, help='Convergence criteria', default=0.0001)
    parser.add_argument('--use_bs', action='store_true', help='set to use bs rather than cphase and camp separately')
    parser.add_argument('--scratch_model', action='store_true', help='set to use a self-generated model rather than FIRST information')
    parser.add_argument('--lotss_file', type=str, help='LoTSS catalogue file', default='' )

    args = parser.parse_args()
    main( args.vis, closure_tels = args.ctels, npix = args.npix, fov_arcsec = args.fov, maxarcmin = args.maxarcmin, zbl = args.zbl, prior_fwhm_arcsec = args.prior_fwhm, doplots = args.doplots, imfile = args.imfile, remove_tels = args.rtels, niter=args.maxiter, cfloor=args.clipfloor, conv_criteria=args.convg, use_bs=args.use_bs, scratch_model=args.scratch_model, do_diagnostic_plots=args.diagplots, lotss_file=args.lotss_file )

