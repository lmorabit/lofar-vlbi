#!/usr/bin/env python
"""
setup for a generic pipeline run of h5_merger.py

written 21 Feb 2022

@author: Leah Morabito)
"""

import os
import glob

def main( h5parm, msin, helperscriptspath_h5merge='', ):

    #msin = msin.split('/')[-1]
    h5parm_out = h5parm.split('/')[-1].replace('.h5','_toapply.h5')

    ## copy the h5_merger.py script
    os.system( 'cp {:s} {:s}'.format( os.path.join( helperscriptspath_h5merge, 'h5_merger.py' ), os.path.join( '.', 'h5_merger.py' ) ) )
    command = '/opt/lofar/pyenv-py2/bin/python h5_merger.py --ms {:s} --h5_tables {:s} --h5_out {:s} --add_cs --h5_time_freq {:s} --circ2lin > convertsols.log'.format( msin, h5parm, h5parm_out, h5parm )
    print(command)
    ## run on the h5parm using a measurement set for the core station info
    os.system( '/opt/lofar/pyenv-py2/bin/python h5_merger.py --ms {:s} --h5_tables {:s} --h5_out {:s} --add_cs --h5_time_freq {:s} --circ2lin > convertsols.log'.format( msin, h5parm, h5parm_out, h5parm ) )

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='cleanup script at the end of the pipeline.')
    parser.add_argument('h5parm')
    parser.add_argument('msin')
    parser.add_argument('--helperscriptspath_h5merge',default='lofar_helpers')

    args = parser.parse_args()

    main(args.h5parm,args.msin,helperscriptspath_h5merge=args.helperscriptspath_h5merge)
