# test script to combine differential and absolute phase tables

# load both tables in
# get a soltab
# use getValues to pull out (phase) values
# combine with numpy
# put back in soltab
# create a new h5parm and write soltab to it

import losoto
from losoto.h5parm import openSoltab, h5parm
from shutil import copyfile
import matplotlib.pyplot as plt
import numpy as np

def plot_phases(h5name):
    phases = openSoltab(h5name, 'sol001', 'phase000').getValues(retAxesVals=False)
    c,_,n,_,p = phases.shape
    plt.plot(phases[0,0,10,0,:])
    plt.show()

def combine_phases(reference,differential,output):

    copyfile(reference,output)
    output = h5parm(output, readonly=False)
    out_soltab = output.getSolset('sol001').getSoltab('phase000')
    diff_soltab = openSoltab(differential, 'sol001', 'phase000')

    ref_phases = out_soltab.getValues(retAxesVals=False)
    print ref_phases.shape
    diff_phases = diff_soltab.getValues(retAxesVals=False)
    
    out_phases = ref_phases + diff_phases 
    # deal with wrapping
    out_phases = np.where(out_phases>np.pi,out_phases-2*np.pi,out_phases)
    out_phases = np.where(out_phases<-np.pi,out_phases+2*np.pi,out_phases)

    out_soltab.setValues(out_phases)
    output.close()
    
if __name__=='__main__':
    plot_phases('sols.h5')
    combine_phases('sols.h5','sols.h5','out.h5')
    plot_phases('out.h5')
    
