#!/usr/bin/env python
import numpy as np
import sys
import os
import multiprocessing
import glob
import logging
import pyrap.tables as pt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse
from astropy.io import ascii
import re

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def closure(vis,tel,lastv=-1,pol=0,use_spw=0,bchan=0,echan=-1,\
            doplot=False,doret=False):

    # Find target source id

    target_id = vis.split('/')[-1].split('_')[0]
    print "target_id", target_id
    print "Antennas for closure phase", tel

    # Find array of requested telescopes and list of telescopes in data

    closure_txt_file = 'closure_text_' + vis
    closure_which_file = 'closure_which_' + vis
    command = 'taql \'select NAME from %s/ANTENNA\' > %s'%(vis,closure_txt_file)
    os.system(command)
    os.system('grep -v select %s >%s'%(closure_txt_file,closure_which_file))
    idxtel = np.loadtxt(closure_which_file,dtype='S')
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
        print 'If no telescopes were found, check that TaQL is installed.'

    aidx_s = np.sort(aidx)

    # Make a smaller MS 'as plain' with the required baseline. This is slow
    # but only needs doing once for an arbitrary number of baselines.
    outfile = 'tmp_'+vis

    command = 'taql \'select from %s where ' % vis
    for i in range (len(aidx_s)):
        for j in range (i+1, len(aidx_s)):
            command += ('ANTENNA1==%d and ANTENNA2==%d' % \
                            (aidx_s[i],aidx_s[j]))
            if i==len(aidx_s)-2 and j==len(aidx_s)-1:
                command += (' giving %s as plain\''% outfile)
            else:
                command += (' or ')

    print 'Selecting smaller MS %s, this will take about 4s/Gb:'% outfile
    os.system (command)

    # Loop around the requested closure triangles

    clstats = np.array([])
    for tr in tel:
        tri = np.array([],dtype='int')
        for i in range(3):
            tri = np.append (tri, aidx[np.argwhere(atel==tr[i])[0][0]])
        tri = np.sort(tri)
        # Make three reference MSs with pointers into the small MS
        command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' %(outfile,tri[0],tri[1],outfile.replace('tmp','tmp_1'))
        os.system(command)
        command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' %(outfile,tri[1],tri[2],outfile.replace('tmp','tmp_2'))
        os.system(command)
        command = 'taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' %(outfile,tri[0],tri[2],outfile.replace('tmp','tmp_3'))
        os.system(command)

        # Load data arrays and get amp, closure phase

        t1 = pt.table(outfile.replace('tmp','tmp_1'))
        t2 = pt.table(outfile.replace('tmp','tmp_2'))
        t3 = pt.table(outfile.replace('tmp','tmp_3'))
        ut = t1.select('TIME')
        spw = t1.select('DATA_DESC_ID')
        d1,d2,d3 = t1.select('DATA'), t2.select('DATA'), t3.select('DATA')
        cp = get_amp_clph(d1[:lastv],d2[:lastv],d3[:lastv],spw[:lastv],
             pol=0, use_spw=use_spw, bchan=bchan, echan=echan)
        try:
            allcp.append(cp)
        except:
            allcp = [cp]

        os.system('rm -fr tmp_*%s'%vis)
	os.system('rm -fr '+closure_txt_file)
	os.system('rm -fr '+closure_which_file)
        clstats = np.append (clstats, np.nanmean(np.gradient(np.unwrap(cp))**2))
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
    if echan==-1 or echan<nchan:
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
    return p1+p2-p3

def addghost (inarray):        # deal with missing frequencies between subbands
    n = len(inarray)
    fdel = []
    for i in range(n):         # first delete if things do not exist
        if not os.path.isdir(inarray[i]):
            print '----->> %s does not exist as an MS'%inarray[i]
            fdel.append(i)
    inarray = np.delete(inarray,fdel)
    for i in range(n):         # make 2-d array of frequencies from chan x files
	spw_table = os.path.join( inarray[i], 'SPECTRAL_WINDOW' )
        newfreq = pt.table(spw_table).getcol('CHAN_FREQ')
        if newfreq.ndim == 2:
            newfreq = newfreq[0]
        try:
            freq = np.vstack((freq,newfreq))
        except:
            freq = np.copy(newfreq)
            fref,chwid = freq[0],freq[1]-freq[0]
            ch_sub = len(freq)

    rfreq=np.rint((freq-fref)/chwid)
    if rfreq.ndim == 1:
        rfreq = np.array([rfreq])    # only one SB, ensure a 2D array
    fo = open('difmap_select','w')   # write select for difmap to exclude flagged if you
    fo.write('select I,')            # run difmap later: difmap crashes on missing chans
    for i in range(len(rfreq)):
        fo.write('%d,%d' % (1+int(rfreq[i,0]),1+int(rfreq[i,-1])))
        fo.write('\n' if i==len(rfreq)-1 else ',')
    fo.close()
    outarray = np.array([inarray[0]])
    for i in range(1,n):        # for each file, see how many ghosts go between
        nghost = int(np.rint((rfreq[i,0]-rfreq[i-1,-1]-1.)/ch_sub))
        for j in range(nghost):
            outarray=np.append(outarray, 'ghost')
        outarray=np.append(outarray,inarray[i])
    return outarray   



def combine_subbands (in1array, nameout, phasecenter, fstep, tstep):
    # get datacolumn
    tmp = phasecenter.split(';')
    datacol = tmp[0]
    phasecenter = tmp[1]
    in2array = addghost(inarray)
    ismissing = False if np.array_equal(in1array,in2array) else True
    fo=open('NDPPP_%s.parset'%nameout,'w')   # write the parset file
    fo.write('msin = [')
    for i in range(len(in2array)):
        fo.write('\'%s\''%in2array[i])
        fo.write(']' if i==len(in2array)-1 else ',')
    fo.write('\n')
    fo.write('msout = '+nameout+'\n')
    fo.write('msout.storagemanager=dysco\n')	
    fo.write('msin.datacolumn = %s\n'%datacol)
    if ismissing:
        fo.write('msin.missingdata=True\n')
        fo.write('msin.orderms=False\n')
    fo.write('steps = [shift,avg,sadder,filter]\n')
    ss = 'shift.phasecenter = [ '+phasecenter.split(',')[0]+'deg,'+phasecenter.split(',')[1]+'deg ]\n'
    fo.write(ss)
    fo.write('shift.type = phaseshift\n')
    fo.write('avg.type = average\n')
    fo.write('avg.timestep = '+str(tstep)+'\n')
    fo.write('avg.freqstep = '+str(fstep)+'\n')
    fo.write('sadder.type = stationadder\n')
    fo.write('sadder.stations = {ST001:\'CS*\'}\n')
    fo.write('filter.type = \'filter\'\n')
    fo.write('filter.baseline = \'!CS*&*\'\n')
    fo.write('filter.remove = True')
    fo.close()
    os.system('NDPPP NDPPP_%s.parset'%nameout)  # run with NDPPP
    mytels = [['DE601','DE605','ST001']]
    scatter_cp = closure( nameout, mytels )
    print '\n Scatter for the direction ' + nameout.replace('.ms','') + ' is %s \n' % scatter_cp
    os.system( 'echo Scatter for the direction '+nameout.replace('.ms','')+' is '+str(scatter_cp)+' >> closure_phases.txt' )
    os.system( 'rm -rf %s'%nameout )
    os.system('rm NDPPP_%s.parset'%nameout)
    

def lotss2coords (lotssfile):

    a = ascii.read(lotssfile)
    if len(a) == 1:
	src = a['Source_id']
	try:
            basestring
	except NameError:
            # Python3, basestring is not longer there.
            basestring = str
	if not isinstance(src, basestring):
            src = 'S' + str(src)
            a = ','.join([str(a['LOTSS_RA']), str(a['LOTSS_DEC']), src])
    coords=np.array([],dtype='S')
    for xx in range(len(a)):
	tmp = a[xx]
	## check if the source_id needs to be converted to a string
	src = tmp['Source_id']
	if type(src) != str:
	    src = 'S'+str(src)
	tmp_1 = ','.join([str(tmp['LOTSS_RA']),str(tmp['LOTSS_DEC']),src])
        coords=np.append(coords,tmp_1)
    return coords


def parallel_function(f,ncpu):            # Scott Sievert parallelize function
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores = int(ncpu)
        print 'Using',ncores,'cores'
        pool = Pool(processes=ncores) # depends on available cores
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallize, f)

def source_thread (i):
    src = i.split(',')[-1]
    print 'PROCESSING SOURCE %s - combining subbands and shifting' % src
    combine_subbands(inarray,src+'.ms',i,8,8)
    print 'PROCESSING SOURCE %s - finished' % src


def source (coords,ncpu):
    source_thread.parallel = parallel_function(source_thread,ncpu)
    parallel_result = source_thread.parallel(coords)

def main( ms_input, lotss_file, ncpu=10, datacol='DATA', nsbs=20 ):

    ncpu = int(ncpu)
    datacol = str(datacol)
    nsbs = int(nsbs)

    if type(ms_input) == str:
        mslist = ms_input.lstrip('[').rstrip(']').replace("'","").replace(" ","").split(',')
    else:
	mslist = ms_input
    mslist = natural_sort( mslist )
    if nsbs < 0:
	nsbs = len(mslist)
    mslist = mslist[0:nsbs]
    print mslist
    global inarray
    inarray = mslist
    coords = lotss2coords( lotss_file )
    tmp = []
    for coord in coords:
	tmp.append( ';'.join([datacol,coord]) )
    coords = tmp
    source( coords, ncpu )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='function to parallel process sources to search for the delay calibrator')
    parser.add_argument('MS_pattern',type=str, help='pattern to search for MS')
    parser.add_argument('lotss_file',type=str,help='catalogue to process')
    parser.add_argument('--ncpu',type=int,help='number of CPUs')
    parser.add_argument('--datacol',type=str,help='datacolumn to use (default DATA)',default='DATA')
    parser.add_argument('--nsbs',type=int,help='number of subbands, default 20, use -1 for all',default=20)
    args = parser.parse_args()

    MS_input = glob.glob( args.MS_pattern )

    main( MS_input, args.lotss_file, ncpu=args.ncpu, datacol=args.datacol, nsbs=args.nsbs )
