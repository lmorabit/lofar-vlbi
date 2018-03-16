#!/usr/bin/env python
import os
import numpy as np
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb as pdb
import copy

###Reading in the the parameters of target data with PYRAP and putting them into directories for further use ###
###Shamelessly borrowed from prefactor pipeline code ###
class ReadMs:
    def __init__(self, ms):
        self.timepara={'start':0, 'end':0, 'step':0, 'cent':0}
	self.freqpara={'start':0, 'end':0, 'step':0, 'cent':0}
	self.msname = ms
	if not os.path.isdir(ms): sys.exit('INPUT MS DOES NOT EXIST!')
        ##########Getting Time parameters first#############
	t = pt.table(ms, readonly=True, ack=False)
	t1 = t.sort ('unique desc TIME')
	self.timepara['step'] = t.getcell('EXPOSURE',0)
	self.timepara['start'] =  np.min(t.getcol('TIME'))-self.timepara['step']/2.
	self.timepara['end'] =  np.max(t.getcol('TIME'))+self.timepara['step']/2.
	self.timepara['cent'] = self.timepara['start']+(self.timepara['end']-self.timepara['start'])/2.
	self.mstimevalues = t1.getcol('TIME')[::-1]
	t1.close()
        ##########Getting Frequency Parameters###################
	freq=pt.table(t.getkeyword("SPECTRAL_WINDOW"))
	self.fullband = freq.getcell('TOTAL_BANDWIDTH', 0)
	self.freqpara['cent'] = freq.getcell('REF_FREQUENCY', 0)
	self.freqpara['step'] = freq.getcell('CHAN_WIDTH', 0)[0]
        self.msfreqvalues = freq.getcell('CHAN_FREQ', 0)
	self.freqpara['start'] = self.msfreqvalues[0]-self.freqpara['step']/2.
	self.freqpara['end'] = self.msfreqvalues[-1]+self.freqpara['step']/2.
	freq.close()	
        ##########Getting Station Names###################
        antennas = pt.table(t.getkeyword("ANTENNA"))
        self.stations = antennas.getcol('NAME')
        antennas.close()
	t.close()
				
    def GetTimepara(self, p=''):
        if p != '': return self.timepara[p]
	else: return self.timepara
    def GetFreqpara(self, p=''):
        if p != '': return self.freqpara[p]
	else: return self.freqpara
    def GetMSNamepara(self): return self.msname
    
def main(parmdbfile, targetms):

    if not os.path.exists(parmdbfile):
        print "Parmdb file %s doesn't exist!" % parmdbfile
        return(1)
    if not os.path.exists(targetms):
        print "Target ms %s doesn't exist!" % targetms
        return(1)
    msinfo = ReadMs(targetms)

    # Open up the parmdb and get station information
    parmdb = pdb.parmdb(parmdbfile)
    parnames = parmdb.getNames()
    parmdb_stations = [ s.split(':')[-1] for s in parnames ]

    ## get an example value for the first CS station
    examplevalue = None
    for name in parnames:
        if "CS" in name:
            csname = name
            break
    csname = csname.split(':')[-1]
    ## get all parameters for the station
    csnames = [ s for s in parnames if csname in s ]
    nparams = len( csnames )

    ## find the first DE station
    for name in parnames:
	if "DE" in name:
	    dename = name
	    break
    dename = dename.split(':')[-1]
    ## get all parameters for the station
    denames = [ s for s in parnames if dename in s ]
    de_nparams = len( denames )

    ## check if the international stations already exist in the right format
    if nparams != de_nparams:   
	## the international stations need to be added
	## get all parameters for the station
	csnames = [ s for s in parnames if csname in s ]
	gain_name = [ s for s in csnames if 'Gain:0:0:Real' in s ][0]
        examplevalue = parmdb.getValuesGrid(gain_name)[gain_name]
	clock_name = [ s for s in csnames if 'Clock:' in s ]
	if len( clock_name ) == 1:
  	    exampleclock = parmdb.getValuesGrid(clock_name[0])[clock_name[0]]
	    exampleclock['values'] = np.zeros(exampleclock['values'].shape)
	else:
	    exampleclock = None
	RotA_name = [ s for s in csnames if 'CommonRotationAngle:' in s ]
	if len( RotA_name ) == 1:
	    exampleRotA = parmdb.getValuesGrid(RotA_name[0])[RotA_name[0]]
	    exampleRotA['values'] = np.zeros(exampleRotA['values'].shape)
	else:
	    exampleRotA = None
	
        # Zero the phases of the example entry
        if examplevalue == None:
            print "Couldn't find an example entry"
            return(1)

        examplevalue['values'] = np.zeros(examplevalue['values'].shape)
        examplevalue_ones = copy.deepcopy(examplevalue)
        examplevalue_ones['values'] = np.ones(examplevalue_ones['values'].shape)

	## find the median value of all core stations
	amps00 = np.zeros(0)
	amps11 = np.zeros(0)
	for name in parnames:
	    if "CS" in name:
		if "0:0:Real" in name:
		    a00 = np.sqrt( parmdb.getValuesGrid(name)[name]['values']**2. + parmdb.getValuesGrid(name.replace('Real','Imag'))[name.replace('Real','Imag')]['values']**2.)
	            if np.count_nonzero(a00) > 0:
	                amps00 = np.append(amps00,a00)
		if "1:1:Real" in name:	
		    a11 = np.sqrt( parmdb.getValuesGrid(name)[name]['values']**2. + parmdb.getValuesGrid(name.replace('Real','Imag'))[name.replace('Real','Imag')]['values']**2.)
	    	    if np.count_nonzero(a11) > 0:
			amps11 = np.append(amps11,a11)

	intl00 = np.median(amps00) * 3.75
	intl11 = np.median(amps11) * 3.75

        # Add the necessary stations
        for antenna_id, antenna in enumerate(msinfo.stations):
            if not "CS" in antenna and not "RS" in antenna:
                ValueHolder = parmdb.makeValue(values=examplevalue['values'],
                                               sfreq=examplevalue['freqs'], 
                                               efreq=examplevalue['freqwidths'],
                                               stime=examplevalue['times'], 
                                               etime=examplevalue['timewidths'], 
                                               asStartEnd=False)
                ValueHolder_ones_00 = parmdb.makeValue(values=examplevalue_ones['values']*intl00,
                                                       sfreq=examplevalue_ones['freqs'], 
                                                       efreq=examplevalue_ones['freqwidths'],
                                                       stime=examplevalue_ones['times'], 
                                                       etime=examplevalue_ones['timewidths'], 
                                                       asStartEnd=False)
                ValueHolder_ones_11 = parmdb.makeValue(values=examplevalue_ones['values']*intl11,
                                                       sfreq=examplevalue_ones['freqs'],
                                                       efreq=examplevalue_ones['freqwidths'],
                                                       stime=examplevalue_ones['times'],
                                                       etime=examplevalue_ones['timewidths'],
                                                       asStartEnd=False)
                ## in case the values already exist (and may be NaN) they need to be deleted first
                parmdb.deleteValues("Gain:0:0:Real:" + antenna)
                parmdb.deleteValues("Gain:1:1:Real:" + antenna)
                parmdb.deleteValues("Gain:0:0:Imag:" + antenna)
                parmdb.deleteValues("Gain:1:1:Imag:" + antenna)
                parmdb.addValues("Gain:0:0:Real:" + antenna,ValueHolder_ones_00)
                parmdb.addValues("Gain:0:0:Imag:" + antenna,ValueHolder)
                parmdb.addValues("Gain:1:1:Real:" + antenna,ValueHolder_ones_11)
                parmdb.addValues("Gain:1:1:Imag:" + antenna,ValueHolder)
		if exampleclock != None:
		    ValueHolder = parmdb.makeValue(values=exampleclock['values'],
                                                   sfreq=exampleclock['freqs'],
                                                   efreq=exampleclock['freqwidths'],
                                                   stime=exampleclock['times'],
                                                   etime=exampleclock['timewidths'],
        	                                   asStartEnd=False)
		    parmdb.deleteValues("Clock:" + antenna)
		    parmdb.addValues("Clock:" + antenna, ValueHolder)
		if exampleRotA != None:
		    ValueHolder = parmdb.makeValue(values=exampleRotA['values'],
                                                   sfreq=exampleRotA['freqs'],
                                                   efreq=exampleRotA['freqwidths'],
                                                   stime=exampleRotA['times'],
                                                   etime=exampleRotA['timewidths'],
                                                   asStartEnd=False)
		    parmdb.deleteValues("CommonRotationAngle:" + antenna)
		    parmdb.addValues("CommonRotationAngle:" + antenna, ValueHolder )


        parmdb.flush()
        parmdb = 0

        return(0)
    else:
	parmdb = 0
	return(0)


if __name__ == "__main__":
    # Check invocation
    print sys.argv[0] + ": modifies a parmdb **in-place** to add international stations with unity gain and zero phase"
    if not (len(sys.argv) == 3 or len(sys.argv) == 4):
        print "Usage: %s <parmdbfile> <targetms>" % sys.argv[0]
        sys.exit()

    # Check that the target files exist
    parmdbfile = sys.argv[1]
    targetms = sys.argv[2]
    main(parmdbfile, targetms)
 
