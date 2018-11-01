import numpy as np,sys,os,multiprocessing,glob
#
# ------- parameters for main script
#
data_dir = '/data020/scratch/lb_bw/long_baseline_pipeline/'
file_prefix = 'L401323_SB'
file_suffix = '_uv.dppp.ndppp_prep_target'
subbands = range(250,350)
ncores = int(os.popen('nproc').read()) - 4
aipsno = 341
indisk = 1
refname = 'ST001'
minmatch = 5
FSTEP,TSTEP = 4,8   # averaging while doing the shifts
ifdivide = 8        # AIPS IFs across band
lotssfile = 'thislotss_lbcs'
flog = 'process.log'
#  -----------------------------------

def imgseq (inna,incl,indisk,inseq):
    uvdata = AIPSUVData (inna,incl,indisk,inseq)
    for i in AIPSCat()[1]:
        if i['name']==inna and i['klass']=='IMGSEQ':
            AIPSImage(inna,'IMGSEQ',indisk,i['seq']).zap()
        if i['name']==inna and i['klass']=='ICL001':
            AIPSImage(inna,'ICL001',indisk,i['seq']).zap()
    imagr=AIPSTask('imagr')
    imagr.indata = uvdata
    imagr.cellsize[1:] = [0.2,0.2]
    imagr.imsize[1:] = [1024,1024]
    imagr.uvwtfn = 'NA'
    imagr.uvrange[1:] = [5.0,100000.0]
    imagr.source[1] = 'BEAM_1'
    imagr.outname = inna
    imagr.niter = 200
    imagr.nchav = 54
    uvtap = [400.0,200.0,100.0,50.0]
    beam  = [0.5,0.9,1.8,5.0]
    for i in range(4):
        imagr.uvtaper[1:] = [uvtap[i],uvtap[i]]
        imagr.bmaj,imagr.bmin = beam[i],beam[i]
        imagr.go()
        AIPSImage(inna,'IBM001',indisk,1).zap()
        AIPSImage(inna,'ICL001',indisk,1).rename(inna,'IMGSEQ',indisk,0)

def lotss2coords (lotssfile):
    a=np.loadtxt(lotssfile,dtype='S')
    if a.ndim==1:
        a=np.array([a])
    coords=np.array([],dtype='S')
    for line in a:
        rahr,ramin,rasec=line[0][4:6],line[0][6:8],line[0][8:12]
        decdeg,decmin,decsec=line[0][13:15],line[0][15:17],line[0][17:21]
        lstr = '%sh%sm%ss,%sd%sm%s'%(rahr,ramin,rasec,decdeg,decmin,decsec)
        coords=np.append(coords,lstr)
    return coords

def addghost (inarray):        # deal with missing frequencies between subbands
    import pyrap; from pyrap import tables as pt
    n = len(inarray)
    fdel = []
    for i in range(n):         # first delete if things do not exist
        if not os.path.isdir(inarray[i]):
            print '----->> %s does not exist as an MS'%inarray[i]
            fdel.append(i)
    inarray = np.delete(inarray,fdel)
    for i in range(n):         # make 2-d array of frequencies from chan x files
        newfreq = pt.table(inarray[i]+'/SPECTRAL_WINDOW').getcol('CHAN_FREQ')
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
    return outarray             # same as input array, but with 'ghost' to make delta_f same

def combine_subbands (in1array, nameout, phasecenter, fstep, tstep):
    inarray = addghost(in1array)
    ismissing = False if np.array_equal(in1array,inarray) else True
    fo=open('NDPPP_%s.parset'%nameout,'w')   # write the parset file
    fo.write('msin = [')
    for i in range(len(inarray)):
        fo.write('\'%s\''%inarray[i])
        fo.write(']' if i==len(inarray)-1 else ',')
    fo.write('\n')
    fo.write('msout = '+nameout+'\n')
    fo.write('msin.datacolumn = DATA\n')
    if ismissing:
        fo.write('msin.missingdata=True\n')
        fo.write('msin.orderms=False\n')
    fo.write('steps = [shift,avg,sadder,filter]\n')
    fo.write('shift.phasecenter = [ '+phasecenter+']\n')
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

def parallel_function(f):            # Scott Sievert parallelize function
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores = max(1,int(os.popen('nproc').read())-5)  # use all cores -5
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
    rn_rh = i.split('h')[0]
    rn_rm = i.split('h')[1].split('m')[0]
    rn_dd = i.split(',')[1].split('d')[0]
    rn_dm = i.split('d')[1].split('m')[0]
    rn = 'L'+rn_rh+rn_rm+'+'+rn_dd+rn_dm
    print 'PROCESSING SOURCE %s - combining subbands and shifting' % rn
    combine_subbands(inarray,rn+'.ms',i,FSTEP,TSTEP)
#   check weights - code snippet from Leah eht_im script
    ss = "taql 'select WEIGHT_SPECTRUM from "+rn+".ms limit 1' > weight_check.txt"
    os.system(ss)
    with open( 'weight_check.txt', 'r' ) as f:
        lines = f.readlines()
    f.close()
    first_weight = np.float(lines[3].lstrip('[').split(',')[0])
    if first_weight < 1e-9:
        ss = "taql 'update "+rn+".ms set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM*1e11'"
        os.system(ss)
        ss = "taql 'update "+rn+".ms set DATA=DATA/1.2e4'"
        os.system(ss)
    print 'PROCESSING SOURCE %s - converting to FITS' % rn
    command = 'ms2uvfits in='+rn+'.ms out='+rn+'.fits writesyscal=F'
    os.system(command)
    print 'PROCESSING SOURCE %s - finished' % rn

def source (coords):
    source_thread.parallel = parallel_function(source_thread)
    parallel_result = source_thread.parallel(coords)
    
#  ----------------------- main script -------------------------

# Part 1 - make and combine measurement sets, write out
#
inarray=[]            # array of measurement sets to read
for i in subbands:
    inarray.append(data_dir+file_prefix+str(i)+file_suffix)

coords = lotss2coords (lotssfile)   # positions to shift to
source(coords)        # make NDPPP shift, average, filter function
fitsfiles = []
for i in coords:
    rn_rh = i.split('h')[0]
    rn_rm = i.split('h')[1].split('m')[0]
    rn_dd = i.split(',')[1].split('d')[0]
    rn_dm = i.split('d')[1].split('m')[0]
    fitsfiles.append('L'+rn_rh+rn_rm+'+'+rn_dd+rn_dm+'.fits')
    
#  Part 2 - AIPS loading and initial flagging, division into IFs for delay cal

from AIPS import AIPS,AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from lofipi_aips import *; from closure_trans import *
from workshop_service import *
AIPS.userno = aipsno
CLIPLEV = 1.2E6

# initial load, clip outliers, form into multiple IFs

for i in fitsfiles:
    inna = i.split('.fits')[0]
    pload('./'+i,inna,indisk,'WPROC')
    pclip(inna,'WPROC',CLIPLEV,indisk)
    puvcop(inna,'WPROC',inna,'UVCOP',1,indisk)
    AIPSUVData(inna,'WPROC',indisk,1).zap()
    pmorif (inna,'UVCOP',inna,'MORIF',ifdivide,indisk)
    AIPSUVData(inna,'UVCOP',indisk,1).zap()

# Part 3 - determine S:N on intl stations from closure phases

antennas = AIPSUVData(inna,'MORIF',1,1).antennas
plsrc = np.array([],dtype='S')
fo = open('clstats.txt','w')
clstats = np.nan*np.ones((len(fitsfiles),len(antennas)))
for i in range(len(antennas)):
    if refname[:minmatch]==antennas[i][:minmatch]:
        refant = i+1
    for j in range(len(fitsfiles)):
        if antennas[i][:2] in ['DE','UK','FR','SE','IE','PL','LT']:
            thissrc = fitsfiles[j].split('.fits')[0]
            if thissrc not in plsrc:
                plsrc = np.append(plsrc,thissrc)
            plotfile = 'cl_'+thissrc+'_'+antennas[i]+'.png'
            triangle = [refname,'RS208',antennas[i]]
            clstats[j,i] = closure_aips(aipsno,triangle,\
                  thissrc,'MORIF',plotfile=plotfile)
            print antennas[i],fitsfiles[j],clstats[j,i]
            fo.write('%s %s %f\n'%( antennas[i],thissrc,clstats[j,i]))

fo.close()

# if we have made some plots, group them

for i in plsrc:
    os.system('montage cl_*%s*.png -geometry 300x300 CL_%s.png'%(i,i))

# need divide with uvmod, multi, indxr here
# Part 4 - full calibration on anything with all telescopes S:N>7

# need to change fullcal for fringe-fitting:
# after fring, run
#   snsmo with doblank 1 and cparm, then snsmo with doblank 0, cparam 0, then
#   sndelay to interpolate over the rest

calsrc = np.array([],dtype='S')
for i in range(len(fitsfiles)):
    if np.nanmin(clstats[i]) > 7:
        thiscal = fitsfiles[i].split('.fits')[0]
        calsrc = np.append (calsrc, thiscal)
        thisra = AIPSUVData(thiscal,'MORIF',indisk,1).header['crval'][3]+360.
        thisdec = AIPSUVData(thiscal,'MORIF',indisk,1).header['crval'][4]
        try:
            calcoord = np.vstack((calcoord,np.array([thisra,thisdec])))
        except:
            calcoord = np.copy(np.array([thisra,thisdec]))
        print '**** DOING FULL CALIBRAION ON %s ****' % thiscal
        fullcal(thiscal, 'MORIF', indisk, refant)

#if calcoord.ndim == 1:
#    calcoord = np.array([calcoord])

# Part 5 - copy the delay + phase solution to the rest of the field

#def skydist (c1,c2):
#    decdist = c1[1]-c2[1]
#    radist = (c1[0]-c2[0])*np.cos(np.deg2rad(c1[1]))
#    return np.hypot(radist,decdist)

#for sfile in fitsfiles:
#    inna = sfile.split('.fits')[0]
#    if inna in calsrc:         # this one is a calibrator, ignore
#        continue
#    coord = np.array([AIPSUVData(inna,'MORIF',indisk,1).table('SU',1)[0]['raobs'],\
#                      AIPSUVData(inna,'MORIF',indisk,1).table('SU',1)[0]['decobs']])
#    imind, mind = -1,1.e9      # find the closest calibrator
#    for i in range(len(calcoord)):
#        thisd = skydist (coord,calcoord[i])
#        if mind > thisd:
#            imind = i
#    ptacop (str(calsrc[imind]),'MORIF',str(inna),'MORIF','CL',0,0,indisk)


        

