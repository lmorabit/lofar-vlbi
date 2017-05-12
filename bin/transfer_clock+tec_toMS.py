#!/usr/bin/env python
import os
import numpy as np
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb as pdb

###Reading in the the parameters of target data with PYRAP and putting them into directories for further use###############
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

# Make an empty parmDB with only the defaults and return the parmdb object
def make_empty_parmdb(outname):
    myParmdb=pdb.parmdb(outname,create=True)
    myParmdb.addDefValues("Gain:0:0:Ampl",1.)
    myParmdb.addDefValues("Gain:1:1:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Ampl",1.)
    myParmdb.addDefValues("Gain:0:0:Real",1.)
    myParmdb.addDefValues("Gain:1:1:Real",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Real",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Real",1.)
    myParmdb.addDefValues("AntennaOrientation",5.497787144)
    myParmdb.addDefValues("RotationMeasure",1e-6)
    return myParmdb

def main(msname, store_basename, npdir='', newparmdbext='-instrument_amp_clock_tec'):

    print 'ENTERED MAIN'
    print msname

    # name (path) for parmdb to be written
    newparmDB = msname+newparmdbext

    # load the numpy arrays written by the previous scripts
    # (filenames constructed in the same way as in these scripts)
#    freqs_ampl = np.load(npdir + '/freqs_for_amplitude_array.npy')
#    amps_array = np.load(npdir + '/'+store_basename + '_amplitude_array.npy')
    clock_array = np.load(npdir + '/fitted_data_dclock_' + store_basename + '_1st.npy')
    dtec_array = np.load(npdir + '/fitted_data_dTEC_' + store_basename + '_1st.npy')
#    freqs_phase = np.load('freqs_for_phase_array.npy')
#    phases_array  = np.load(store_basename + '_phase_array.npy')
#    station_names = np.load(npdir + '/'+store_basename + '_station_names.npy')

    #print "phases shape:",np.shape(phases_array)
    #print "amps shape:",np.shape(amps_array)
    #print "clock shape:",np.shape(clock_array)

    #for ms in mslist: #this script works only on one MS!
    msinfo = ReadMs(msname)

    station_names = msinfo.stations
    ## get a list of core stations
    cs = [ sn for sn in station_names if 'CS' in sn ]
    ## get rid of core stations and replace with ST001
    station_names = [ sn for sn in station_names if 'CS' not in sn ]
    station_names.append('ST001')


    # this is the same for all antennas
    starttime = msinfo.timepara['start']
    endtime   = msinfo.timepara['end']
    startfreqs = msinfo.msfreqvalues-msinfo.GetFreqpara('step')/2.
    endfreqs   = msinfo.msfreqvalues+msinfo.GetFreqpara('step')/2.
    ntimes  = 1
    nfreqs  = len(startfreqs)

    outDB = make_empty_parmdb(newparmDB)

    # Now do the interpolating
    for antenna_id, antenna in enumerate(station_names):
        if antenna not in msinfo.stations:
            pass

        #handle the clock-value (no fancy interpolating needed)
        clock_pdb = np.array( np.median(clock_array[:,antenna_id]) ,ndmin=2)
        ValueHolder = outDB.makeValue(values=clock_pdb,
                                      sfreq=startfreqs[0], efreq=endfreqs[-1],
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Clock:'+antenna,ValueHolder)

	#handle the TEC-value (no fancy interpolating needed)
        dtec_pdb = np.array( np.median(dtec_array[:,antenna_id]), ndmin=2)
	ValueHolder = outDB.makeValue(values=dtec_pdb,
				      sfreq=startfreqs[0], efreq=endfreqs[-1],
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('TEC:'+antenna,ValueHolder)


        ## add core stations if it's the ST001 station
        if antenna == 'ST001':
	    ## clock values
            ValueHolder = outDB.makeValue(values=clock_pdb,
                                          sfreq=startfreqs[0], efreq=endfreqs[-1],
                                          stime=starttime, etime=endtime, asStartEnd=True)
	    for cs_ant in cs:
		outDB.addValues('Clock:'+cs_ant,ValueHolder)
	
	    ## tec values
	    ValueHolder = outDB.makeValue(values=dtec_pdb,
                                          sfreq=startfreqs[0], efreq=endfreqs[-1],
                                          stime=starttime, etime=endtime, asStartEnd=True)
            for cs_ant in cs:
                outDB.addValues('TEC:'+cs_ant,ValueHolder)

    outDB = False
    return {'transfer_parmDB': newparmDB }
