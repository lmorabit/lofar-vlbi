#!/usr/bin/env python
import sys,os,time,numpy as np,pyrap,matplotlib
from pyrap import tables as pt
from matplotlib import pyplot as plt

#  closure.py v1  Neal Jackson 2017 Jan 4
# Given a visibility file, return the scatter on the closure phases for a
# high signal-to-noise triangle. Scatter is 1.64 for random closure phases
# and about 1.00 for reasonably strong sources.
# modified to work with the genericpipeline version of the Long Baseline pipeline
# Alexander Drabent, 2017 Jan 12
# modified to select only 2 MHz of bandwidth if the file is larger
#
# Neal Jackson, 2018 Mar 19  v3
# modified so that the input station array is a list of 3xn stations, separated
# by semicolons. Each set of 3 will be extracted from the MS and the program
# returns an array of closure statistics. (If only one triangle, a float
# is returned in order to be backwards compatible with version 2). Optionally
# the closure phases are plotted, or the closure phases themselves are returned.
# For efficiency reasons, first cuts out the required baselines from the
# large data file into a separate measurement set and only then indexes
# individual baselines into the smaller dataset. (I have done some timing
# and think that this is the fastest way to do this).
#
# NJ 2018 May 19 v4
#
#   modified. instead of producing closure phase at the input time resolution,
#   will average by factors of up to 30 and produce the best closure phase
#   statistic. This is the best obtainable statistic with data averaging,
#   limited by the time at which the coherence on the individual baselines
#   is bad. See Lofar-LB memo number 3.
#
# arguments:
#    vis = visibility file name
#    tel = np.array nx3 strings of telescopes to use
#    lastv = only process this number of data points (default all)
#    pol = polarization (default first in file, RR)
#    use_spw = 0 SPW to use, default 0
#    bchan, echan: start and end channels within this spectral window
#    doplot: if True produce plot of the closure phases
#    doret: if True return a list of closure phase arrays in the order
#           requested in the tel input
#
#  returns:
#  Most general case, dopipe=False,doret=True: will return a list, one 
#  entry for each MS where each entry is itself a list of closure phase 
#  scatter, followed (if doret) by numpy arrays of closure phases

def closure(vis,tel,lastv=-1,pol=0,use_spw=0,bchan=0,echan=-1,\
            doplot=False,doret=False):

    # Find target source id

    target_id = vis.split('/')[-1].split('_')[0]
    print "target_id", target_id
    print "Antennas for closure phase", tel

    # Find array of requested telescopes and list of telescopes in data

    command = 'taql \'select NAME from %s/ANTENNA\' >closure_txt'%vis
    os.system(command)
    os.system('grep -v select closure_txt >closure_which')
    idxtel = np.loadtxt('closure_which',dtype='S')
    atel = np.unique(np.ravel(tel))

    # For each requested telescope, determine its position in the list
    # If more than one telescope in the list match to within the number
    #   of letters in the requested telescope, keep the first (so CS002
    #   will match to CS002HBA0 and CS002HBA1 will be ignored)
    # Keep a list of telescopes not found, to print if we need to crash
    
    notfound = []
    aidx = np.array([],dtype='int')
    for a in atel:
        found_this = False
        for i in range(len(idxtel)):
            if a==idxtel[i][:len(a)]:
                aidx = np.append(aidx,i)
                found_this = True
        if not found_this:
            notfound.append (a)

    if len(notfound):
        print 'The following telescopes were not found:',notfound


    aidx_s = np.sort(aidx)

    # Make a smaller MS 'as plain' with the required baseline. This is slow
    # but only needs doing once for an arbitrary number of baselines.

    if os.path.exists ('cl_temp.ms'):
        os.system('rm -fr cl_temp.ms')
    command = 'taql \'select from %s where ' % vis
    for i in range (len(aidx_s)):
        for j in range (i+1, len(aidx_s)):
            command += ('ANTENNA1==%d and ANTENNA2==%d' % \
                            (aidx_s[i],aidx_s[j]))
            if i==len(aidx_s)-2 and j==len(aidx_s)-1:
                command += (' giving cl_temp.ms as plain\'')
            else:
                command += (' or ')

    print 'Selecting smaller MS cl_temp.ms, this will take about 4s/Gb:'
    os.system (command)

    # Loop around the requested closure triangles

    clstats = np.array([])
    for tr in tel:
        tri = np.array([],dtype='int')
        for i in range(3):
            tri = np.append (tri, aidx[np.argwhere(atel==tr[i])[0][0]])
        tri = np.sort(tri)
        # Make three reference MSs with pointers into the small MS
        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp1.ms\'' %(tri[0],tri[1])
        os.system(command)
        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp2.ms\'' %(tri[1],tri[2])
        os.system(command)
        command = 'taql \'select from cl_temp.ms where ANTENNA1==%d and ANTENNA2==%d giving closure_temp3.ms\'' %(tri[0],tri[2])
        os.system(command)

        # Load data arrays and get amp, closure phase

        t1 = pt.table('closure_temp1.ms')
        t2 = pt.table('closure_temp2.ms')
        t3 = pt.table('closure_temp3.ms')
        ut = t1.select('TIME')
        spw = t1.select('DATA_DESC_ID')
        d1,d2,d3 = t1.select('DATA'), t2.select('DATA'), t3.select('DATA')
        clthis, cp = get_amp_clph(d1[:lastv],d2[:lastv],d3[:lastv],spw[:lastv],
             pol=0, use_spw=use_spw, bchan=bchan, echan=echan)
        try:
            allcp.append(cp)
        except:
            allcp = [cp]

        os.system('rm -fr closure_temp*ms')
        if os.path.exists ('closure_which'):
            os.system('rm closure_which')
        clstats = np.append (clstats, clthis)
    if doplot:
        cl_mkplot (allcp,tel,target_id)
    clstats = clstats[0] if len(clstats)==1 else clstats
    if doret:
        return clstats,allcp
    else:
        return clstats

def cl_mkplot(allcp,tel,target_id):
    ny = int(np.floor(np.sqrt(np.float(len(allcp)))))
    nx = 1+len(allcp)/ny if len(allcp)%ny else len(allcp)/ny
    matplotlib.rcParams.update({'font.size':8-nx//2})
    for i in range(len(allcp)):
        ax=plt.subplot(nx,ny,i+1)
        plt.plot(allcp[i],'b,')
        plt.xlabel('Sample number')
        plt.ylabel('Phase/rad')
        plt.text(0.01,0.9,'%s %s-%s-%s'%(target_id,tel[i,0],tel[i,1],tel[i,2]),\
                 transform=ax.transAxes)
    plt.savefig('%s_closure.png'%target_id,bbox_inches='tight')

def get_amp_clph(d1,d2,d3,spw,pol=0,use_spw=0,bchan=0,echan=-1):
    a1,a2,a3,cp = np.array([]),np.array([]),np.array([]),np.array([])
    p1,p2,p3 = np.array([]),np.array([]),np.array([])
    nchan = np.int(len(d1[0]['DATA']))
    bchan = max(bchan,0)
    if echan==-1 or echan>nchan:  # bugfix
        echan=nchan
    for i in range(len(d1)):
        if spw[i].values()[0] != use_spw:
            continue
        vis1,vis2,vis3 = d1[i]['DATA'][bchan:echan,pol],\
                         d2[i]['DATA'][bchan:echan,pol],\
                         d3[i]['DATA'][bchan:echan,pol]
        pd1,pd2,pd3 = np.nansum(vis1),np.nansum(vis2),np.nansum(vis3)
        p1 = np.append(p1, np.arctan2 (pd1.imag,pd1.real))
        p2 = np.append(p2, np.arctan2 (pd2.imag,pd2.real))
        p3 = np.append(p3, np.arctan2 (pd3.imag,pd3.real))
    clph = p1+p2-p3
    clph = clph[~np.isnan(clph)]
    np.putmask (clph,clph<-np.pi,clph+2*np.pi)
    np.putmask (clph,clph>np.pi,clph-2*np.pi)
    x,y,z = np.arange(1,30),np.array([]),np.array([])
    for i in x:
        y = np.append(y,np.mean(np.gradient(np.unwrap(phavg(clph,i)))**2))
    pfit = np.polyfit (x,y,3)
    for i in x:
        z = np.append(z,np.poly1d(pfit)(i))
    return z.min(),clph

def phavg (phase, n):
    phase = phase[0:n*(len(phase)//n)]
    rreal,rimag = np.cos(phase), np.sin(phase)
    arreal = np.average (np.reshape(rreal,(-1,n)),axis=1)
    arimag = np.average (np.reshape(rimag,(-1,n)),axis=1)
    return np.arctan2 (arimag, arreal)

#ctel = 'ST001;DE601;DE603;ST001;FR606;DE609'
#closure ('L328370.ms',ctel,doplot=True)


def main(ms_input,station_input,lastv=-1,pol=0,use_spw=0,bchan=0,\
             echan=-1,doplot=False,doret=False,dopipe=True):
    """
    Deriving closure phases of all directions
   
    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    station_input : List of stations for closure triangles, separated by
        semicolons. There must be nx3 of those.
    lastv (int):   Index of last visibility
    pol     (int, default 0): Polarization number
    use_spw (int, default 0): Spectral window number
    bchan   (int, default 0): First channel number
    echan   (int, default -1): Last channel number
    doplot  (bool, default False): If true, generate plot
    doret   (bool, default False): If true, return list of closure phase
               arrays in second argument (DO NOT USE IF dopipe=TRUE)
    dopipe  (bool, default True): For backwards compatibility. If dopipe
            is True, will write text in the closure_phase file as before.
            If True, station_input will be truncated to 3 stations.
        
    Returns
    -------
    result : dict
        Dict with the name of the generated mapfile
    """
    
    mslist = str(ms_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(',')
    stationlist = str(station_input).lstrip('[').rstrip(']').replace("'","").replace(" ","").split(';')
    tel = np.asarray(stationlist)
    if len(tel)%3:
        print 'Station list is not a multiple of 3 long, truncating it...'
        tel = tel[0:3*(len(tel)//3)]
    if dopipe and len(tel)>3:
        print 'dopipe True, truncating to 3 telescopes'
        tel = tel[:3]
    tel = tel.reshape(len(tel)/3,3)
    for ms in mslist:
        print 'Now operating on', ms
        scatter_cp = closure(ms,tel,lastv=lastv,pol=pol,\
                             use_spw=use_spw,bchan=bchan,echan=echan,\
                             doplot=doplot,doret=doret)
        if dopipe:
            print '\n Scatter for the direction ' + ms.split('/')[-1].split('_')[0] + ' is %s \n' % scatter_cp
            os.system('echo Scatter for the direction ' + ms.split('/')[-1].split('_')[0] + ' is ' + str(scatter_cp) + ' >> closure_phases.txt')
        else:
            try:
                allret.append(scatter_cp)
            except:
                allret = [scatter_cp]

    if not dopipe:
        return allret
