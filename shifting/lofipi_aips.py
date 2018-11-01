from math import *
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from scipy import ndimage; from scipy.ndimage import measurements
import re,sys,pickle,numpy as np,os,glob,time,warnings; from numpy import fft

def ptacop (name1, class1, name2, class2, ttype, innum, outnum, disk):
    tacop = AIPSTask('TACOP')
    tacop.indata = AIPSUVData(name1, class1, disk, 1)
    tacop.outdata = AIPSUVData(name2, class2, disk, 1)
    tacop.inext = ttype
    tacop.invers = innum
    tacop.outvers = outnum
    tacop.go()
    
def pwtmod (aipsname, refant, antennas, supweight=50.0, indisk=1):
    uvdata = AIPSUVData (aipsname, 'FITS', indisk, 1)
    wtmod = AIPSTask ('WTMOD')
    wtmod.indata = uvdata
    wtmod.outdata = uvdata
    wtmod.antwt[1:] = [1]*len(antennas)
    wtmod.antwt[antennas.index(refant)+1] = supweight
    wtmod.inp()
    wtmod.go()

# dparm(8): binary with 1=rates, 2=delays, 4=phase

def pfring (aipsname,refant,antennas,source,indisk=1,delaywin=600,ratewin=20,\
            solint=1,snr=2,logdir='./',weightit=3,zero=0,aipsclass='FITS',\
            aipsseq=1,suppress_rate=0,in2name='',in2class='',in2seq=1,\
            docalib=-1):
    uvdata = AIPSUVData (aipsname,aipsclass,indisk,aipsseq)
    fring = AIPSTask ('FRING')
    fring.refant = refant
    fring.indata = uvdata
    fring.calsour[1:] = [source]
    fring.antennas[1:] = antennas
    if in2name!='' and in2class=='':
        uv2data = AIPSUVData (in2name, in2class, indisk, in2seq)
        fring.in2data = uv2data
    fring.solint = solint
    fring.aparm[1:] = [0,0,0,0,0,2,snr,0,0]
    fring.dparm[1:] = [0,delaywin,ratewin,0,0,0,0,zero,suppress_rate]
    fring.weightit = weightit
    fring.docalib = docalib
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
#    fring.inp()
    fring.go()
    sys.stdout.close(); sys.stdout = stdout

def pclcal (aipsname,indisk,inver,aipsclass='FITS',logdir='./',\
            snver=-1,gainver=0,gainuse=0):
    uvdata = AIPSUVData (aipsname,aipsclass,indisk,1)
    clcal = AIPSTask ('clcal')
    clcal.indata = uvdata
    clcal.inver = inver
    clcal.snver = inver if snver==-1 else snver
    clcal.gainver = gainver
    clcal.gainuse = gainuse
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    clcal.go()
    sys.stdout.close(); sys.stdout = stdout

def psplit (aipsname,source,indisk,logdir='./'):
    uvdata = AIPSUVData (aipsname,'FITS',indisk,1)
    split = AIPSTask ('split')
    split.indata = uvdata
    split.outclass = 'SPLIT'
    split.docalib = 1
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    split.go()
    sys.stdout.close(); sys.stdout = stdout
    uvdata = AIPSUVData(source,'SPLIT',indisk,1)
    uvdata.rename(aipsname,'SPLIT',1)

def pload (filename,aipsname,indisk,outcl,logdir='./',doindxr=True,\
           idxint=5./60.):
    fitld = AIPSTask ('FITLD')
    fitld.datain = str(filename)
    fitld.outna = aipsname
    fitld.outcl = outcl
    fitld.dokeep = 1
    fitld.outdisk = indisk
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    fitld.go ()
    if doindxr:
        uvdata = AIPSUVData (aipsname,outcl,1,1)
        cltables = np.array([],dtype='int')   # remove extant CL tables
        for i in range(len(uvdata.tables)):
            if uvdata.tables[i][1] == 'AIPS CL':
                np.append (cltables,uvdata.tables[i][0])
        for i in cltables[::-1]:
            print 'Removing existing CL table %d'%i
            uvdata.table('CL',i).zap()
        indxr = AIPSTask ('INDXR')     # INDXR creating new table
        indxr.cparm[1:] = [0,0,float(idxint),0,0,0,0,0,0,0]
        indxr.indata = uvdata
        indxr.go()
    sys.stdout.close(); sys.stdout = stdout

def stars (aipsname, incl, indisk, intext='./starsfile',logdir='./'):
    stars = AIPSTask('stars')
    stars.inname = aipsname
    stars.inclass = incl
    stars.indisk = indisk
    try:
        stars.stvers = 0    # does not exist in some AIPS versions
    except:
        pass
    stars.intext = './starsfile'
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
#    stars.inp()
    try:
        stars.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def greys (aipsname, incl, indisk, pmin, pmax, stfac, stvers, logdir='./'):
    greys = AIPSTask('greys')
    greys.inname = aipsname
    greys.inclass = incl
    greys.indisk = indisk
    greys.pixrange[1:] = [float(pmin),float(pmax)]
    greys.dotv = -1
    greys.dowedge = -1
    greys.ltype = 7
    greys.stfac = stfac
    try:
        greys.stvers = stvers  # does not exist in some aips versions
    except:
        pass
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    try:
        greys.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def lwpla (aipsname,incl,indisk,outfile,logdir='./'):
    lwpla = AIPSTask('lwpla')
    lwpla.inname = aipsname
    lwpla.inclass = incl
    lwpla.indisk = indisk
    lwpla.outfile = outfile
    lwpla.docolor = 1
    lwpla.plcolors[10][1] = 1
    lwpla.plcolors[10][2] = 1
    lwpla.plcolors[10][3] = 1
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    try:
        lwpla.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def pimagr (aipsname,aipsclass,docalib,source='',imsize=256,cellsize=0.1,\
            nchav=-1,disk=1,gainuse=0,robust=0,niter=100,outname='',\
            antennas = [], stokes='I',uvwtfn=''):
    imagr = AIPSTask('imagr')
    imagr.inname = aipsname
    imagr.inclass = aipsclass
    imagr.docalib = docalib
    imagr.indisk = disk
    imagr.stokes = stokes
    imagr.imsize[1:] = [imsize,imsize]
    imagr.cellsize[1:] = [cellsize,cellsize]
    if len(antennas):
        imagr.antennas[1:] = antennas
#        imagr.baseline[1:] = antennas
    if source=='':
        try:
            su = AIPSUVData(aipsname,aipsclass,disk,1).table('SU',1)
            imagr.source[1] = su[0]['source']
            print 'Extracted source %s' % imagr.source[1]
        except:
            print 'Unable to extract source from source table'
    if nchav==-1:
        try:
            h = AIPSUVData(aipsname,aipsclass,disk,1).header
            imagr.nchav = h['naxis'][h['ctype'].index('FREQ')]
        except:
            imagr.nchav=0
    imagr.gainuse = gainuse
    imagr.robust = robust
    imagr.niter = niter
    imagr.outname = outname
#    imagr.inp()
    imagr.go()

def pclip(aipsname,aipsclass,cliplevel,indisk):
    clip = AIPSTask('clip')
    clip.inname = aipsname
    clip.inclass = aipsclass
    clip.indisk = indisk
    clip.aparm[1] = float(cliplevel)
    clip.aparm[2] = float(cliplevel)
    clip.aparm[3] = 1.0e-6
    clip.aparm[4] = 1.0e-6
    clip.go()

def puvcop(aipsname,aipsclass,outname,outclass,flagver,indisk):
    uvcop = AIPSTask('uvcop')
    uvcop.inname = aipsname
    uvcop.inclass = aipsclass
    uvcop.indisk = indisk
    uvcop.outname = outname
    uvcop.outclass = outclass
    uvcop.outdisk = indisk
    uvcop.flagver = flagver
    uvcop.go()

def pmorif(aipsname,aipsclass,outname,outclass,npiece,indisk):
    morif = AIPSTask('morif')
    morif.inname = aipsname
    morif.inclass = aipsclass
    morif.indisk = indisk
    morif.outname = outname
    morif.outclass = outclass
    morif.outdisk = indisk
    morif.npiece = npiece
    morif.go()

def psnsmo(aipsname,aipsclass,indisk,smotype,boxcar=2,delaycut=50,smooth=1):
    snsmo = AIPSTask('snsmo')    
    snsmo.inname = aipsname
    snsmo.inclass = aipsclass
    snsmo.indisk = indisk
    snsmo.smotype = smotype
    if smotype=='DELA' or smotype=='VLBI':
        snsmo.bparm[4] = snsmo.bparm[5] = snsmo.bparm[9] = snsmo.bparm[10] = smooth
        snsmo.cparm[4] = snsmo.cparm[5] = boxcar
        snsmo.cparm[9] = snsmo.cparm[10] = delaycut
        snsmo.doblank = -1
    snsmo.go()

def pcalib(aipsname,aipsclass,indisk,solint,refant,docalib=-1,solmode='P',soltype='L1',calsour='BEAM_1',weightit=1,in2name='',in2class='',in2seq=0,invers=0,ncomp=-10000,snr=3):
    calib = AIPSTask('calib')
    calib.calsour[1] = calsour
    calib.inname = aipsname
    calib.inclass = aipsclass
    calib.indisk = indisk
    calib.in2name = in2name
    calib.in2class = in2class
    calib.in2seq = in2seq
    calib.docalib = docalib
    calib.solmode = solmode  # e.g. 'P', 'A&P'
    calib.soltype = soltype
    calib.solint = solint
    calib.weightit = weightit
    calib.aparm[7] = snr
    calib.go()

    
    
