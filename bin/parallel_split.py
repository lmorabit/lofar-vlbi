#!/usr/bin/env python
import numpy as np
import sys
import os
import multiprocessing
import glob
import logging
import casacore.tables as casatb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse
from astropy.table import Table
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
	## explicitly open and close the table
	with casatb.table(spw_table,readonly='True') as t:
	    newfreq = t.getcol('CHAN_FREQ')
	t.close()
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

def combine_subbands (inarray, nameout, datacol, phasecenter, fstep, tstep, phscmd, filcmd, nthreads):

    # get datacolumn
    in2array = addghost(inarray)
    if np.array_equal(inarray,in2array):
	ismissing = False
    else:
	ismissing = True
    parset_name = 'NDPPP_%s.parset'%nameout.replace('MS','').replace('ms','')
    with open(parset_name,'w') as fo:
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
        fo.write('numthreads=%s\n'%str(nthreads))
#        fo.write('steps = [shift,avg,sadder,filter]\n')
        fo.write('steps = [shift,avg]\n')
        fo.write('shift.phasecenter = [%s]\n'%phasecenter)
        fo.write('shift.type = phaseshift\n')
        fo.write('avg.type = average\n')
        fo.write('avg.timestep = '+str(tstep)+'\n')
#        fo.write('avg.freqstep = '+str(fstep)+'\n')
#        fo.write('sadder.type = stationadder\n')
#        fo.write("sadder.stations = %s\n"%phscmd)
#        fo.write('filter.type = \'filter\'\n')
#        fo.write("filter.baseline = %s\n"%filcmd)
#        fo.write('filter.remove = True')
    fo.close()
    os.system('NDPPP %s'%parset_name)  # run with NDPPP
    os.system('rm NDPPP_%s'%parset_name) # remove the parset
    
def lotss2coords (lotssfile):

    a = Table.read(lotssfile)
    coords=np.array([],dtype='S')
    for xx in range(len(a)):
	tmp = a[xx]
	## check if the source_id needs to be converted to a string
	src = tmp['Source_id']
	if type(src) != str:
            if str(src)[0:1] == 'I':
                src = str(src)
            else:
                src = 'S'+str(src)
	tmp_1 = ','.join([str(tmp['RA']),str(tmp['DEC']),src])
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

    ## unpack the dictionary
    src = i['src']
    print 'PROCESSING SOURCE %s - combining subbands and shifting' % src
    inarray = i['inarray']
    datacol = i['datacol']
    freqstep = i['freqstep']
    timestep = i['timestep']
    phasecen = i['phasecen']    
    phaseupcmd = i['phaseup_cmd']
    filtercmd = i['filter_cmd']
    nameout = i['outname']
    nthreads = i['nthreads']
    combine_subbands(inarray,nameout,datacol,phasecen,freqstep,timestep,phaseupcmd,filtercmd,nthreads)
    print 'PROCESSING SOURCE %s - finished' % src


def source (coords,ncpu):
    source_thread.parallel = parallel_function(source_thread,ncpu)
    parallel_result = source_thread.parallel(coords)

def main( ms_input, lotss_file, phaseup_cmd="{ST001:'CS*'}", filter_cmd='!CS*&&*', ncpu=10, datacol='DATA', timestep=8, freqstep=8, nsbs=999, nthreads=0 ):

    phaseup_cmd = str(phaseup_cmd)
    filter_cmd = str(filter_cmd)
    ncpu = int(ncpu)
    nthreads = int(nthreads)
    if nthreads == 0:
	nthreads = int(4)
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
    print( 'check 2' )
    ## get the calibrator names
    cal_names = []
    for coord in coords:
	cal_names.append(coord.split(',')[2])

    ## check if you need to concat bands first
    if nsbs < len(mslist):
	## need to split into bands
	print('YOU HAVE CHOSEN POORLY')
    else:
	print('COMBINING ALL SUBBANDS')
	# define an empty dictionary
	cal_dict = {}
	# add the mslist
	cal_dict['inarray'] = mslist
        ## get obsid
        tmp = mslist[0].split('/')[-1]
        obsid = tmp.split('_')[0]
	cal_dict['obsid'] = obsid
        coords = lotss2coords( lotss_file )
	## add information like datacolumn, etc.
	cal_dict['datacol'] = datacol
	cal_dict['freqstep'] = str(freqstep)
	cal_dict['timestep'] = str(timestep)
	cal_dict['phaseup_cmd'] = phaseup_cmd
	cal_dict['filter_cmd'] = filter_cmd
        cal_dict['nthreads'] = nthreads

        tmp = []
        for coord in coords:
	    # create a temporary dictionary
	    tmp_dict = cal_dict.copy()
	    coord_tmp = coord.split(',')
	    tmp_dict['src'] = coord_tmp[2]
            phasecen = ','.join([coord_tmp[0]+'deg',coord_tmp[1]+'deg'])
	    tmp_dict['phasecen'] = phasecen
	    tmp_dict['outname'] = coord_tmp[2]+'_'+obsid+'_imdir.ms'
	    tmp.append(tmp_dict)
        coords = tmp
	## find number of sources to set right number of cpus
	## each ndppp process will use 4 threads (this is hardcoded atm)
        if len(coords) > int(ncpu/nthreads):
	    ncpu = int( len(coords) )
        else:
	    ncpu = int( ncpu/nthreads )
        print( 'HELLOOOOOOOO' )
        print( 'using ncpu: ', ncpu )

        #print( coords )
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
    parser.add_argument('--filter_cmd',type=str,default='!CS*&&*')
    parser.add_argument('--ncpu',type=int,help='number of CPUs',required=True)
    parser.add_argument('--nthreads',type=int,help='number of threads',required=True,default=4)
    parser.add_argument('--datacol',type=str,help='datacolumn to use (default CORRECTED_DATA)',default='CORRECTED_DATA')
    parser.add_argument('--timestep',type=int,default=8)
    parser.add_argument('--freqstep',type=int,default=8)
    parser.add_argument('--nsbs',type=int,help='number of subbands to combine before combining all (default 999=all)',default=999)
    args = parser.parse_args()

    MS_input = glob.glob( args.MS_pattern )

    main( MS_input, args.lotss_file, phaseup_cmd=args.phaseup_cmd, filter_cmd=args.filter_cmd, ncpu=args.ncpu, datacol=args.datacol, timestep=args.timestep, freqstep=args.freqstep, nsbs=args.nsbs )
