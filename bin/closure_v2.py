#!/usr/bin/env python
import sys
import os
import numpy as np
import pyrap; from pyrap import tables as pt
import matplotlib
from matplotlib import pyplot as plt
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

#  closure.py v1  Neal Jackson 2017 Jan 4
# Given a visibility file, return the scatter on the closure phases for a
# high signal-to-noise triangle. Scatter is 1.64 for random closure phases
# and about 1.00 for reasonably strong sources.
# modified to work with the genericpipeline version of the Long Baseline pipeline
# Alexander Drabent, 2017 Jan 12
# modified to select only 2 MHz of bandwidth if the file is larger

def closure(vis, tel, lastv=-1):
    # Find target source id
    target_id = vis.split('/')[-1].split('_')[0]
    print "target_id", target_id
    print "Antennas for closure phase", tel
    # Find which number is which antenna (amend if you are a taql superstar and can get this less clanky)
    command = 'taql \'select NAME from %s/ANTENNA\' >closure_which'%vis
    os.system(command)
    idx_tels, iline = [-1,-1,-1], 0
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            for i in range(3):
                if tel[i]==line[:len(tel[i])]:
                    idx_tels[i] = iline
            iline += 1
    f.close()
    if -1 in idx_tels:
        print 'Did not find one or more of the telescopes'
        exit()

    # get the bandwidth
    command = 'taql \'select TOTAL_BANDWIDTH from %s/SPECTRAL_WINDOW\' >total_bandwidth'%vis
    os.system(command)
    with open('total_bandwidth') as f:
	lines = f.readlines()
    f.close()
    os.system('rm total_bandwidth')
    bandwidth = np.float(lines[-1].rstrip('\n'))
    if bandwidth > 2e6:
	# get channel width
        command = 'taql \'select CHAN_WIDTH from %s/SPECTRAL_WINDOW\' >chan_width'%vis
	os.system(command)
        with open('chan_width') as f:
            lines = f.readlines()
        f.close()
        os.system('rm chan_width')
        chan_width = np.float(lines[-1].split(',')[0].strip('['))
	if chan_width < 195312:
	    nchans = np.ceiling(195312/chan_width)
        else:
	    nchans = 1


    # Make three reference MSs with pointers
    os.system('rm -fr closure_temp*.ms')
    command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving closure_temp1.ms\'' %(vis,idx_tels[0],idx_tels[1])
    os.system(command)
    command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving closure_temp2.ms\'' %(vis,idx_tels[1],idx_tels[2])
    os.system(command)
    command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving closure_temp3.ms\'' %(vis,idx_tels[0],idx_tels[2])
    os.system(command)
    # Load data arrays (NB assumes all three have the same structure - needs changing)
    t1 = pt.table('closure_temp1.ms')
    t2 = pt.table('closure_temp2.ms')
    t3 = pt.table('closure_temp3.ms')
    ut = t1.select('TIME')
    spw = t1.select('DATA_DESC_ID')
    d1,d2,d3 = t1.select('DATA'), t2.select('DATA'), t3.select('DATA')
    cp,a1,a2,a3 = get_amp_clph(d1[:lastv],d2[:lastv],d3[:lastv],spw[:lastv])
    # mask wrapped values
    # unwrap phases: it will correct an array of phases modulo 2pi such that all jumps are less than or equal to pi
    cp = np.unwrap(cp)
    # return a statistic - is 1.64 for random closure phase, less for coherent
    os.system('rm -fr closure_temp*ms')
    os.system('rm closure_which')
    os.system('rm closure_out_%s.png' % target_id)
    plt.clf()
    plt.ylabel(r'$\phi$ $[rad]$') ; plt.plot(cp,'k+') #; plt.grid(True)
    plt.title('Scatter for target %s is %s' % (target_id,np.mean(np.gradient(cp)**2)))
    plt.savefig('closure_out_%s.png' % target_id)
    return np.mean(np.gradient(cp)**2)

def get_amp_clph(d1,d2,d3,spw,nspw=0,pol=0,nchans=1):
    a1,a2,a3,cp = np.array([]),np.array([]),np.array([]),np.array([])
    p1,p2,p3 = np.array([]),np.array([]),np.array([])
    for i in range(len(d1)):
        if spw[i].values()[0] != nspw:
            continue
        vis1,vis2,vis3 = d1[i]['DATA'][0:1,pol],d2[i]['DATA'][0:1,pol],d3[i]['DATA'][0:1,pol]
        a1 = np.append(a1,abs(vis1).sum()/len(vis1))
        a2 = np.append(a2,abs(vis2).sum()/len(vis2))
        a3 = np.append(a3,abs(vis3).sum()/len(vis3))
        p1 = np.append(p1,np.arctan2 (vis1.sum().imag,vis1.sum().real))
        p2 = np.append(p2,np.arctan2 (vis2.sum().imag,vis2.sum().real))
        p3 = np.append(p3,np.arctan2 (vis3.sum().imag,vis3.sum().real))
    return p1+p2-p3,a1,a2,a3


def main(ms_input, station_input):
    """
    Deriving closure phases of all directions
    
    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
        
    Returns
    -------
    result : dict
        Dict with the name of the generated mapfile
    """
    
    mslist = str(ms_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(',')
    stationlist = str(station_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(';')
    os.system('rm closure_phases.txt')
    os.system('touch closure_phases.txt')
    for ms in mslist:
      print 'Now operating on', ms
      scatterclosurephase = closure(ms, stationlist, lastv=-1)
      print '\n Scatter for the direction ' + ms.split('/')[-1].split('_')[0] + ' is %s \n' % scatterclosurephase
      os.system('echo Scatter for the direction ' + ms.split('/')[-1].split('_')[0] + ' is ' + str(scatterclosurephase) + ' >> closure_phases.txt')
      pass
    
  
    pass
