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
import datetime

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

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
    outarray = np.array([inarray[0]])
    for i in range(1,n):        # for each file, see how many ghosts go between
        nghost = int(np.rint((rfreq[i,0]-rfreq[i-1,-1]-1.)/ch_sub))
        for j in range(nghost):
            outarray=np.append(outarray, 'ghost')
        outarray=np.append(outarray,inarray[i])
    return outarray   

def combine_subbands (in1array, nameout, datacol, phasecenter, fstep, tstep, phscmd, filcmd):

    # get datacolumn
    in2array = addghost(inarray)
    ismissing = False if np.array_equal(in1array,in2array) else True
    fo=open('NDPPP_%s'%nameout.replace('MS','parset'),'w')   # write the parset file
    fo.write('msin = [')
    for i in range(len(in2array)):
        fo.write('\'%s\''%in2array[i])
        fo.write(']' if i==len(in2array)-1 else ',')
    fo.write('\n')
    fo.write('msout = '+nameout+'\n')
    fo.write('msin.datacolumn = %s\n'%datacol)
    if ismissing:
        fo.write('msin.missingdata=True\n')
        fo.write('msin.orderms=False\n')
    fo.write('numthreads=3\n')
    fo.write('steps = [shift,avg,sadder,filter]\n')
    ss = 'shift.phasecenter = [%s]\n'%phasecenter
    fo.write(ss)
    fo.write('shift.type = phaseshift\n')
    fo.write('avg.type = average\n')
    fo.write('avg.timestep = '+str(tstep)+'\n')
    fo.write('avg.freqstep = '+str(fstep)+'\n')
    fo.write('sadder.type = stationadder\n')
    fo.write("sadder.stations = %s\n"%phscmd)
    fo.write('filter.type = \'filter\'\n')
    fo.write("filter.baseline = %s\n"%filcmd)
    fo.write('filter.remove = True')
    fo.close()
    os.system('NDPPP NDPPP_%s'%nameout.replace('MS','parset'))  # run with NDPPP
    os.system('rm NDPPP_%s'%nameout.replace('MS','parset')) # remove the parset
    
def lotss2coords (lotssfile):

    a = ascii.read(lotssfile)
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
    tmp = i.split(';')
    src = tmp[-1]
    print 'PROCESSING SOURCE %s - combining subbands and shifting' % src
    datacol = tmp[0]
    freqstep = tmp[1]
    timestep = tmp[2]
    phasecen = tmp[3]
    phaseupcmd = tmp[4]
    filtercmd = tmp[5]
    combine_subbands(inarray,src,datacol,phasecen,freqstep,timestep,phaseupcmd,filtercmd)
    print 'PROCESSING SOURCE %s - finished' % src


def source (coords,ncpu):
    source_thread.parallel = parallel_function(source_thread,ncpu)
    parallel_result = source_thread.parallel(coords)

def main( ms_input, lotss_file, phaseup_cmd="{ST001:'CS*'}", filter_cmd='!CS*&*', ncpu=10, datacol='DATA', timestep=8, freqstep=8, nsbs=10 ):

    phaseup_cmd = str(phaseup_cmd)
    filter_cmd = str(filter_cmd)
    ncpu = int(ncpu)
    datacol = str(datacol)
    timestep = int(timestep)
    freqstep = int(freqstep)
    nsbs = int(nsbs)

    if type(ms_input) == str:
        mslist = ms_input.lstrip('[').rstrip(']').replace("'","").replace(" ","").split(',')
    else:
	mslist = ms_input
    mslist = natural_sort( mslist )

    ## housekeeping
    coords = lotss2coords( lotss_file )
    ## get the calibrator names
    cal_names = []
    for coord in coords:
	cal_names.append(coord.split(',')[2])

    ## check if you need to concat bands first
    if nsbs < len(mslist):
	## need to split into bands
	nbands = np.ceil(np.float(len(mslist))/nsbs)
	for xx in np.arange(nbands):
	    if len(mslist) >= nsbs:
                bandlist = mslist[0:(nsbs)-1]
	        mslist = mslist[nsbs:-1]
	    else:
	    	bandlist = mslist

   	    global inarray
            inarray = bandlist
	    ## get obsid
	    tmp = bandlist[0].split('/')[-1]
            obsid = tmp.split('_')[0]
            coords = lotss2coords( lotss_file )
            tmp = []
            for coord in coords:
                coord_tmp = coord.split(',')
                phasecen = ','.join([coord_tmp[0]+'deg',coord_tmp[1]+'deg'])
                tmp.append(';'.join([datacol,str(freqstep),str(timestep),phasecen,phaseup_cmd,filter_cmd,coord_tmp[2]+'_'+obsid+'_phasecal_band'+str(xx)+'.MS']))
            coords = tmp
            starttime = datetime.datetime.now()
            source( coords, ncpu )
            endtime = datetime.datetime.now()
            timediff = endtime - starttime
            print( 'total time (sec): %s'%str( timediff.total_seconds() ) )
 
	## now combine the bands
        for cal_name in cal_names:
	    cal_bands = glob.glob(cal_name+'*phasecal_band*MS')
	    cal_bands = natural_sort(cal_bands)
	    cal_bands = addghost( cal_bands )
	    tmp = cal_bands[0].split('_')
	    nameout = '_'.join([tmp[0],tmp[1],tmp[2]])
	    nameout = nameout + '.MS'
	    parset_name = 'NDPPP_%s_total.parset'%cal_name
	    with open(parset_name,'w') as fo:  # write the parset file
                fo.write('msin = [')
                ss = ""
                for cal_band in cal_bands:
                    ss = ss + "'%s',"%cal_band
                ss = ss.rstrip(',')
                ss = ss + ']'
                fo.write(ss+'\n')
                fo.write('msout = '+nameout+'\n')
		fo.write('numthreads=3\n')
                fo.write('steps = []\n')
	    fo.close()
    	    os.system('NDPPP %s'%parset_name)  # run with NDPPP
            os.system('rm %s'%parset_name) # remove the parset
	    ## remove the intermediate bands
	    for cal_band in cal_bands:
	        os.system('rm -r %s'%cal_band )

    else:
        global inarray
        inarray = mslist
        ## get obsid
        tmp = mslist[0].split('/')[-1]
        obsid = tmp.split('_')[0]
        coords = lotss2coords( lotss_file )
        tmp = []
        for coord in coords:
	    coord_tmp = coord.split(',')
            phasecen = ','.join([coord_tmp[0]+'deg',coord_tmp[1]+'deg'])
	    tmp.append(';'.join([datacol,str(freqstep),str(timestep),phasecen,phaseup_cmd,filter_cmd,coord_tmp[2]+'_'+obsid+'_phasecal.MS']))
        coords = tmp
        print( coords )
        starttime = datetime.datetime.now()
        source( coords, ncpu )
        endtime = datetime.datetime.now()
        timediff = endtime - starttime
        print( 'total time (sec): %s'%str( timediff.total_seconds() ) )



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='function to parallel process sources to search for the delay calibrator')
    parser.add_argument('MS_pattern',type=str, help='pattern to search for MS')
    parser.add_argument('lotss_file',type=str,help='catalogue to process')
    parser.add_argument('--phaseup_cmd',type=str,default="{ST001:'CS*'}")
    parser.add_argument('--filter_cmd',type=str,default='!CS*&*')
    parser.add_argument('--ncpu',type=int,help='number of CPUs',required=True)
    parser.add_argument('--datacol',type=str,help='datacolumn to use (default CORRECTED_DATA)',default='CORRECTED_DATA')
    parser.add_argument('--timestep',type=int,default=8)
    parser.add_argument('--freqstep',type=int,default=8)
    parser.add_argument('--nsbs',type=int,help='number of subbands to combine before combining all (default 10)',default=10)
    args = parser.parse_args()

    MS_input = glob.glob( args.MS_pattern )

    main( MS_input, args.lotss_file, phaseup_cmd=args.phaseup_cmd, filter_cmd=args.filter_cmd, ncpu=args.ncpu, datacol=args.datacol, timestep=args.timestep, freqstep=args.freqstep, nsbs=args.nsbs )
