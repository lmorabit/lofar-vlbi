#!/usr/bin/env python
import os
import glob
import argparse

def main( ms_input, h5pattern='S*_tec.h5' ):

    tec_h5parms = glob.glob(h5pattern)

    for tec_h5parm in tec_h5parms:

        filestem = tec_h5parm.replace('.h5','')
        losoto_parset = filestem + '_losoto.parset' 
        try:
            imstem = filestem.split('_')[0] + '_' + filestem.split('_')[1] + '_' + filestem.split('_')[2] + '_'
        except:
            imstem = filestem
        with open( losoto_parset, 'w' ) as f:
            f.write('[plotTEC]\n')
            f.write('plotFlag = False\n')
            f.write('axesInPlot = [time]\n')
            f.write('prefix = %s\n'%imstem)
            f.write('soltab = sol000/tec000\n')
            f.write('axisInTable = ant\n')
            f.write('operation = PLOT\n')
            f.write('refAnt = ST001')
        f.close()

        os.system( 'losoto -V %s %s > %s 2>&1'%(tec_h5parm, losoto_parset, losoto_parset.replace('parset','log') ) )
        ## remove parset and log
	os.system( 'rm {:s}'.format(losoto_parset) )
	os.system( 'rm {:s}'.format(losoto_parset.replace('parset','log') ) )



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='plot tec solutions for all existing h5parms')
    parser.add_argument('MS_pattern',type=str, help='pattern to search for MS')
    parser.add_argument('--h5_pattern',type=str, default='S*_tec.h5', help='pattern to search for h5parms')
    args = parser.parse_args()

    main( args.MS_pattern, h5pattern=args.h5_pattern )
