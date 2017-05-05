#!/usr/bin/env python

# converted from locapiv3b.py by Colm Coughlan, March 2016

import os
import sys

# short locapi script to run mscorpol on current file
# (mscorpol doesn't have a main function suitable for direct pipeline integration)
# can also use other convertors

def main(msname, cpath, mode='lin2circ'):

	if mode == 'lin2circ':
		os.system('python '+cpath+' -i '+msname + ' -c DATA -o DATA -p True')
	else:
		if mode == 'mscorpol':
			os.system('python '+path+' -f '+msname)	# execute mscorpol command
		else:
			if mode == 'simple2circ':
				os.system('python '+cpath+' '+msname)
			else:
				raise ValueError('Error, unknown mode: '+mode)

    
if __name__ == "__main__":
    # Options
    print('Mostly a pythonplugin, but attempting a run anyway...')
    print('Input form: msname, path to converter, mode: lin2circ, mscorpol, simple2circ')
    if(len(sys.argv)==4):
		main(sys.argv[1], sys.argv[2],sys.argv[3])
    else:
		print('Error: check your inputs.')


