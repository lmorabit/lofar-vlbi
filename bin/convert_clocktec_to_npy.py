#!/usr/bin/env python
import sys
import os
import numpy
import pyrap.tables
from pylab import *
import sys, os, glob, re
import numpy as np
import shutil
#import progressbar
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm
#import lofar.expion.fitting as fitting

args = sys.argv
globaldbname=args[1]
calsource=args[2]
########################
###### USER INPUT ######

#globaldbname = 'L128487.h5' # input h5 parm file
#calsource    = '3C295' # name for writing outputfiles

#### END USER INPUT ####
########################


pi = numpy.pi
c  = 2.99792458e8
ionmodel = h5parm(globaldbname ,readonly=True)
ionsolset = ionmodel.getSolset('sol000')

## get the individual tables
amptab = ionsolset.getSoltab('amplitude000')
phasetab = ionsolset.getSoltab('phase000')
phaseOffsetTab = ionsolset.getSoltab('phase_offset000')
clocktab = ionsolset.getSoltab('clock000')
tectab = ionsolset.getSoltab('tec000')
anttab = ionsolset.getAnt()
print('tables read in ...')

source_id     = 0  # source ID in global_db (usually 0)
refantenna_id = 0

## frequencies for phase array
freqs = np.copy(phasetab.freq)
subbands = np.unique(np.round(freqs/195.3125e3-512.))
nsubbands = len(subbands)
nchan = len(freqs)/nsubbands
if nsubbands*nchan != len(freqs):
    print "convert_clocktec_to_npy.py: iregular number of ch/SB detected! Bailing out!"
    print "  nchan %d, nSB: %d, nfreq: %d" % (nchan, nsubbands, len(freqs))
    sys.exit(1)
tmpfreqs = freqs.reshape([nsubbands,nchan])
freq_per_sb = np.mean(tmpfreqs,axis=1)

print('... converting tables to npy ...')
amparray = np.asarray( amptab.val )
phasearray = np.asarray( phasetab.val )
## the offset is constant in time and frequency
phaseOffsetArray = np.asarray( phaseOffsetTab.val )
## extend default values to other frequencies
tmpPhaseOffset = phaseOffsetArray.repeat( len(freq_per_sb) )
tmpPhaseOffset2 = tmpPhaseOffset.reshape( len(anttab), len(freq_per_sb) )
phaseOffsetArray = np.transpose( tmpPhaseOffset2 )
clockarray = np.asarray( clocktab.val )
tecarray   = np.asarray( tectab.val )


## station names 
stationsnames = [ stat for stat in phasetab.ant]

print('... writing files ...')
## clock and tec values
os.system('rm -f ' + 'fitted_data_dclock_' + calsource + '_1st.sm.npy')
os.system('rm -f ' + 'fitted_data_dTEC_'   + calsource + '_1st.sm.npy')
#numpy.save('fitted_data_dclock_' + calsource + '_1st.sm.npy', clockarray)
numpy.save('fitted_data_dclock_' + calsource + '_1st.sm.npy', clockarray)
#numpy.save('fitted_data_dTEC_'   + calsource + '_1st.sm.npy', tecarray)
numpy.save('fitted_data_dTEC_'   + calsource + '_1st.sm.npy', tecarray)

## phase information
os.system('rm -f freqs_for_phase_array.npy')
os.system('rm -f ' + calsource + '_phase_array.npy')
os.system('rm -f ' + calsource + '_station_names.npy')
np.save('freqs_for_phase_array.npy', freq_per_sb)
np.save(calsource + '_phase_array.npy', phaseOffsetArray)
np.save(calsource + '_station_names.npy', stationsnames)
print('... files written.')

