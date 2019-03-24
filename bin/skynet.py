#!/usr/bin/env python
import numpy as np
import os
import sys
import glob
from astropy.coordinates import SkyCoord
from astropy.io import ascii


def write_skymodel (ra,dec,model,outname):

    print model

    if outname!='':
        f = open(outname,'w')
        f.write ('# (Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, Orientation) = format\n')
    for i in range(len(model)):
	if isinstance( ra, str ):
	    sra = ra
	    sdec = dec
	else:
	    cosd = 3600.*np.cos(np.deg2rad(dec))
            s = SkyCoord(ra-model[i,0]/cosd,dec+model[i,1]/3600,unit='degree')
            s = s.to_string(style='hmsdms')
            sra = s.split()[0]
	    sdec = s.split()[1]
        sra = sra.replace('h',':').replace('m',':').replace('s','')
        sdec = sdec.replace('d','.').replace('m','.').replace('s','')
        print 'ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f'%(i,sra,sdec,model[i,2],\
              model[i,3],model[i,3]*model[i,4],np.rad2deg(model[i,5]))
        if outname!='':
            f.write('ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f\n'%(i,sra,sdec,model[i,2],\
                  model[i,3],model[i,3]*model[i,4],np.rad2deg(model[i,5])))
    if outname!='':
        f.close()

################## skynet ##############################

def main (vis, delayCalFile='' ):

    ## make sure the parameters are the correct format
    vis = vis.rstrip('/')
    vis = vis.split('/')[-1]
    vis_src = vis.split('_')[0]

    ## get flux from primary_delay_calibrator.csv
    a = ascii.read(delayCalFile)
    for xx in range(len(a)):
        tmp = a[xx]
        src = tmp['Source_id']
        if type(src) != str:
            src = 'S'+str(src)
	if src == vis_src:
	    ra = tmp['LOTSS_RA']
	    dec = tmp['LOTSS_DEC']
	    smodel = tmp['Total_flux']

    print 'generating point model'
    point_model = np.array( [ [0.0,0.0,smodel,0.1,0.0,0.0] ] )
    write_skymodel (ra,dec,point_model,vis+'/skymodel')
    # run makesourcedb to generate sky
    ss = 'rm -fr %s/sky\n'%vis
    print ss
    os.system( ss )
    ss = 'makesourcedb in=%s/skymodel out=%s/sky format=\'<\''%(vis,vis)
    print ss
    os.system ( ss )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Skynet script to handle LBCS calibrators.")

    parser.add_argument('vis', type=str, help='Measurement set for which to run skynet')
    parser.add_argument('--delay_cal_file',type=str,help='delay calibrator information')

    args = parser.parse_args()

    main( args.vis, delayCalFile=args.delay_cal_file )

