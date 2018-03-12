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
    
def main(parmdbfile, targetms, phaseonly = True):

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
    ## check if the international stations exist already
    station_check = [ s for s in msinfo.stations if s in parmdb_stations ]
    if len( station_check ) < len( msinfo.stations ):   
	# if there are fewer stations in the parmdb than the measurement set
        examplevalue = None
        for name in parnames:
            if "CS" in name:
                examplevalue = parmdb.getValuesGrid(name)[name]
                break

        # Zero the phases of the example entry
        if examplevalue == None:
            print "Couldn't find an example entry"
            return(1)

        examplevalue['values'] = np.zeros(examplevalue['values'].shape)
        if not(phaseonly):
            examplevalue_ones = copy.deepcopy(examplevalue)
            examplevalue_ones['values'] = np.ones(examplevalue_ones['values'].shape)

	## find the median value of all core stations
	amps00 = np.zeros(0)
	amps11 = np.zeros(0)
	for name in parnames:
	    if "CS" in name:
		if "0:0:Real" in name:
		    a00 = np.sqrt( parmdb.getValuesGrid(name)[name]['values']**2. + parmdb.getValuesGrid(name.replace('Real','Imag'))[name.replace('Real','Imag')]['values']**2.)
		if "1:1:Real" in name:	
		    a11 = np.sqrt( parmdb.getValuesGrid(name)[name]['values']**2. + parmdb.getValuesGrid(name.replace('Real','Imag'))[name.replace('Real','Imag')]['values']**2.)
	    if np.count_nonzero(a00) > 0:
		amps00 = np.append(amps00,a00)
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

                if phaseonly:
    		    ## in case the values already exist (and may be NaN) they need to be deleted first
                    parmdb.deleteValues("Gain:0:0:Phase:" + antenna)
                    parmdb.deleteValues("Gain:1:1:Phase:" + antenna)	
                    parmdb.addValues("Gain:0:0:Phase:" + antenna,ValueHolder)
                    parmdb.addValues("Gain:1:1:Phase:" + antenna,ValueHolder)
                else:
                    ValueHolder_ones = parmdb.makeValue(values=examplevalue_ones['values'],
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
                    parmdb.addValues("Gain:0:0:Real:" + antenna,ValueHolder_ones*intl00)
                    parmdb.addValues("Gain:0:0:Imag:" + antenna,ValueHolder)
                    parmdb.addValues("Gain:1:1:Real:" + antenna,ValueHolder_ones*intl11)
                    parmdb.addValues("Gain:1:1:Imag:" + antenna,ValueHolder)

        parmdb.flush()
        parmdb = 0

        return(0)
   else:
	parmdb = 0
	return(0)


if __name__ == "__main__":
    # Check invocation
    print sys.argv[0] + ": modifies a phase-only or amp and phase parmdb **in-place** to add international stations with unity gain and zero phase"
    if not (len(sys.argv) == 3 or len(sys.argv) == 4):
        print "Usage: %s <parmdbfile> <targetms> <optional:1 for phaseonly parmdb, 0 for amp and phase>" % sys.argv[0]
        sys.exit()

    # Check that the target files exist
    parmdbfile = sys.argv[1]
    targetms = sys.argv[2]
    if(len(sys.argv) == 4):
        do_phase = int(sys.argv[3])
        if do_phase == 1:
            do_phase = True
        else:
            do_phase = False
    else:
        do_phase = True
    main(parmdbfile, targetms, do_phase)
 
