#!/usr/bin/env python
"""
setup for a generic pipeline run of facetselfcal.py

written 21 Feb 2022

@author: Leah Morabito)
"""

import os
import glob

def main( msin, helperscriptspath='', helperscriptspath_h5merge='', configfile='', destdir='' ):

    msin = msin.split('/')[-1]
    skymod = os.path.join( msin, 'skymodel' )

    ## make a new directory and copy the measurement set into it
    os.mkdir('delay_solve')
    destdir = os.path.join( destdir, 'delay_solve' )
    os.system( 'cp -r {:s} delay_solve/'.format( msin ) )
    os.chdir('delay_solve')

    with open( 'run_selfcal.log', 'a+') as f:
        f.write('{:s}\n{:s}\n{:s}\n{:s}\n{:s}\n{:s}'.format(msin,skymod,helperscriptspath,helperscriptspath_h5merge,configfile,destdir) )

    # copy the config file
    destfile = os.path.join( destdir, configfile.split('/')[-1] )
    ss = 'cp {:s} {:s}'.format(configfile, destfile)
    print( ss )
    os.system( ss )

    # update it
    os.system( 'sed -i "s~FACETSELFCAL_DIR~{:s}~g" {:s}'.format(helperscriptspath,destfile) )
    os.system( 'sed -i "s~LOFARHELPERS_DIR~{:s}~g" {:s}'.format(helperscriptspath_h5merge,destfile) )
    os.system( 'sed -i "s~MYMODEL~{:s}~g" {:s}'.format(skymod,destfile) )

    ## also copy the h5_merger.py script
    os.system( 'cp {:s} {:s}'.format( os.path.join( helperscriptspath_h5merge, 'h5_merger.py' ), os.path.join( destdir, 'h5_merger.py' ) ) )

    os.system( 'cp {:s} {:s}'.format( os.path.join( helperscriptspath, 'polconv.py' ), os.path.join( destdir, 'polconv.py' ) ) )

    os.system( '/opt/lofar/pyenv-py2/bin/python {:s} {:s}'.format(os.path.join(helperscriptspath,'facetselfcal.py'), msin ) )

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='cleanup script at the end of the pipeline.')
    parser.add_argument('msin')
    parser.add_argument('--helperscriptspath',default='lofar_facet_selfcal')
    parser.add_argument('--helperscriptspath_h5merge',default='lofar_helpers')
    parser.add_argument('--configfile',default='facetselfcal_config.txt')
    parser.add_argument('--destdir',default='jobdir')

    args = parser.parse_args()

    main(args.msin,helperscriptspath=args.helperscriptspath,helperscriptspath_h5merge=args.helperscriptspath-h5merge,configfile=args.configfile,destdir=args.destdir)
