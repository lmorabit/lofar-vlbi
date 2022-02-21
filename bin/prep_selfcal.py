#!/usr/bin/env python
# -* coding: utf-8 -*-

"""
setup for a generic pipeline run of facetselfcal.py

written 21 Feb 2022

@author: Leah Morabito)
"""

import argparse
import os
import glob
import logging

def main(msin,skymod,,helperscriptspath='',helperscriptspath-h5merge='',configfile='',destdir=''):

    # copy the config file
    destfile = os.path.join( destdir, configfile.split('/')[-1] )
    ss = 'cp {:s} {:s}'.format(configfile, destfile)
    print( ss )
    os.system( ss )

    # update it
    os.system( 'sed -i "s~FACETSELFCAL_DIR~{:s}~g" {:s}'.format(helperscriptspath,destfile) )
    os.system( 'sed -i "s~LOFARHELPERS_DIR~{:s}~g" {:s}'.format(helperscriptspath-h5merge,destfile) )
    os.system( 'sed -i "s~MYMODEL~{:s}~g" {:s}'.format(skymod,destfile) )
    with open( destfile, 'r' ) as f:
        lines = f.readlines()
    print( lines )

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='cleanup script at the end of the pipeline.')
    parser.add_argument('msin')
    parser.add_argument('skymod')
    parser.add_argument('--helperscriptspath',default='lofar_facet_selfcal')
    parser.add_argument('--helperscriptspath-h5merge',default='lofar_helpers')
    parser.add_argument('--configfile',default='facetselfcal_config.txt')
    parser.add_argument('--destdir',default='jobdir')

    args = parser.parse_args()

    main(args.msin,args.skymod,helperscriptspath=args.helperscriptspath,helperscriptspath-h5merge=args.helperscriptspath-h5merge,configfile=args.configfile,destdir=args.destdir)
