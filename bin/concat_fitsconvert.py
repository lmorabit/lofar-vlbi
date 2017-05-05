#!/usr/bin/env python
# -*- coding: utf-8 -*-

# converted from locapiv3b.py by Colm Coughlan, March 2016

import os, sys
import pyrap.tables as pt

# locapi script to concat together grouped subbands into final dataset. Convert concatenated file to fits

def main(file_list, ms2uvfits_path='ms2uvfits'):
    file_list = file_list.lstrip('[').rstrip(']').split(',')
    print file_list
    outstem = (file_list[0].split('.')[0])+'.ms'
    
    pt.msconcat (file_list, outstem)	# virtually concat files
    os.system(ms2uvfits_path+' in='+outstem+' out='+outstem.replace('.ms','.fits')+' writesyscal=F') 

    return {'fits': outstem.replace('.ms','.fits') }	# return name of final file created to main pipeline for further use
    
if __name__=='__main__':
    # Options
    print('Mostly a pythonplugin, but attempting a run anyway...')
    print('Input form: input mapfile, path to ms2uvfits (optional)')
    if(len(sys.argv)==3):		
		main(sys.argv[1], sys.argv[2])
    else:
        if(len(sys.argv)==2):		
		    main(sys.argv[1])
        else:
		    print('Error: check your inputs.')


