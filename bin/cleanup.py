#!/usr/bin/env python
# -* coding: utf-8 -*-

"""
Cleanup temporary files after running long baseline pipeline

Created on Tue Aug 28 2018

@author: Leah Morabito)
"""

import argparse
import os
import glob
import logging

def main():

    # delete prep target files
    os.system('rm S*ndppp_prep_target')
    # and the individual h5parms
    os.system('rm S*.h5')

    # get a list of concatenated files
    concat_files = glob.glob('S*msdpppconcat')
    lengths = []
    for concat_file in concat_files:
	lengths.append(len(concat_file))
    # delete them if they aren't the longest file name
    for xx in range(len(lengths)):
	if xx < max(lengths):
	    os.system('rm -r %s'%concat_files[xx] )

    # rename delay calibrator something sensible
    delay_cal_file = glob.glob('S*msdpppconcat')[0]
    tmp = delay_cal_file.split('_')
    new_filename = tmp[0] + '_' + tmp[1] + '_delaycal.msdpppconcat'
    os.system( 'mv %s %s'%(delay_cal_file, new_filename) )


if __name__ == "__main__":

    main()
