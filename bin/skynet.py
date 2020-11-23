#!/usr/bin/env python
import numpy as np
import os
import sys
import glob
from astropy.coordinates import SkyCoord
from astropy.table import Table


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
    tmp = vis.split('/')[-1]
    vis_src = tmp.split('_')[0]

    ## get flux from best_delay_calibrators.csv
    t = Table.read( delayCalFile, format='csv' )
    ## find the RA column
    mycols = t.colnames
    ra_col = [ val for val in mycols if val == 'RA' ]
    de_col = [ val for val in mycols if val == 'DEC' ]
    if len(ra_col) == 1:
        ra_col = ra_col[0]
        de_col = de_col[0]
    else:
        ## pick LoTSS position
        ra_col = [ val for val in mycols if val == 'RA_LOTSS' ][0]
        de_col = [ val for val in mycols if val == 'DEC_LOTSS' ][0]
    ## get the right source
    src_ids = t['Source_id']

    src_names = []
    for src_id in src_ids:
	if type( src_id ) != str:
            if str(src_id)[0:1] == 'I':
		val = str(src_id)
            else:
		if vis_src[0:1] == 'S':
                    val = 'S'+str(src_id)
		else:
		    val = str(src_id)
	else:
	    val = src_id
        src_names.append(val)
    src_idx = [ i for i, val in enumerate(src_names) if vis_src == val ][0]

    ra = t[ra_col].data[src_idx]
    dec = t[de_col].data[src_idx]
    smodel = t['Total_flux'].data[src_idx]*1.0e-3

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

