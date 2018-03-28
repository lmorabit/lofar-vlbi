#!/usr/bin/env python
import numpy as np
import os
import sys
import glob
import astropy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
import time
import pyrap.tables as pt
# neal specific
import scipy
from scipy import fftpack as ff
import warnings
import multiprocessing
from scipy import ndimage,optimize
from astropy.coordinates import SkyCoord


# Requires:
# 1) Measurement set
# 2) list of three closure telescopes (e.g. ['TS001','DE601','DE605'])

# For the measuement set, evaluates the closure-scatter. If this is <CTHR, makes a model

# modes: 
# 1) make a short-baseline image (phases already calibrated) and selfcal to that
# 2) start with a point source
# 3) start with something from the model-engine

# Returns: closure threshold value
#          Images will be left in directory depending on CTh and mode

# 4ch/8s 20Gb per source (FOV) 1ch/16s 3 GB/source to go to 5' fields
# eor scripts gives out cal table as a lofar parmdb

###################### taql_funcs ################################
def taql_calc (vis, vistable, qtab, qtype):
    os.system('taql \'CALC '+qtype+' ([select '+qtab+' from '+vis+\
              '/'+vistable+'])\' >taql_out')
    f=open('taql_out')
    for v in f:
        try:
            val = float(v.rstrip('\n'))
        except:
            pass
    f.close()
    os.system('rm taql_out')
    return val

def taql_num (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f=open('taql_out')
    for v in f:
        if 'select result of' in v:
            n = int(v.split('of')[1].split('row')[0])
            break
    f.close()
    os.system('rm taql_out')
    return n

def taql_from (vis,vistable,qtab):
    os.system('taql \'select '+qtab+' from '+vis+'/'+vistable+'\' >taql_out')
    f = open('taql_out')
    for v in f:
        pass
    f.close()
    os.system('rm taql_out')
    return v.rstrip('\n').rstrip(']').lstrip('[').split(',')

################## idx_tels #############################


# From a MS, return the index numbers of telescopes. Neal Jackson 03/17
# Input telescopes given as an array of strings
# Note that the index numbers are found in one of two places, the NAME and
# STATION columns (latter if it has been written by AIPS, for instance), so
# this routine deals with that too. NB returns the OFFSET (starts at 0) not
# the AIPS telescope number (starts at 1) in this case.
# Output is a list of indices of telescopes. Note that you need to use the
# lower one first if reading data out of a MS.

def get_idx_tels (data, tel):
    os.system('taql \'select NAME from %s/ANTENNA\' >closure_which'%data)
    f = open('closure_which')
    for line in f:
        if 'select' in line or not len(line.rstrip('\n')):
            continue
        try:
            a = int(line) # it's just a number, use STATION column instead
            f.close()
            os.system('taql \'select STATION from %s/ANTENNA\' >closure_which'%data)
            break
        except:
            f.close()
            break
    idx_tels, iline = [-1]*len(tel), 0
    f = open('closure_which')
    for line in f:
        if not 'select' in line:
            for i in range(len(tel)):
                if tel[i]==line[:len(tel[i])]:
                    idx_tels[i] = iline
            iline += 1
    f.close()
    if -1 in idx_tels:
        os.system('cat closure_which')
        print 'Did not find one or more of the telescopes'
        print 'Telescopes present are those in list above'
        return []
    os.system('rm closure_which')
    return idx_tels

################## mkgauss ##############################

def mkgauss (naxes,pos,flux,fwhm,axrat=1.0,angle=0.0,ignore=4.0,dodist=False):
# note that total flux = peak flux in a pixel * 1.1331*FWHM**2
# angle is major axis East of North
    a = np.zeros (naxes[0]*naxes[1]).reshape(naxes[1],naxes[0])
    fwhm /= 1.66667
    if axrat==1.0 and angle==0.0:
        for i in range (naxes[1]):
            ydist=float(i)-pos[1]
            for j in range (naxes[0]):
                xdist=float(j)-pos[0]
                if xdist*xdist+ydist*ydist>ignore*ignore*fwhm*fwhm:
                    continue
                if not dodist:
                    a[i,j] = flux*np.exp(-(xdist*xdist+ydist*ydist)/ \
                                    (fwhm*fwhm))/(fwhm*fwhm*np.pi)
                else:
                    a[i,j] = np.hypot(xdist,ydist)
        return a
    sinth = np.sin(angle*np.pi/180.0)
    costh = np.cos(angle*np.pi/180.0)
    r = np.array([-sinth,costh,-costh,-sinth])
    rt = np.array([-sinth,-costh,costh,-sinth])
    sig = np.array([fwhm,0.0,0.0,fwhm*axrat])
    scr1 = mxmul (sig,r)
    scr2 = mxmul (rt, scr1)
    scr1 = mxinv (scr2)
    for i in range(naxes[1]):
        ydist=float(i)-pos[1]
        if abs(ydist)>ignore*fwhm:
            continue
        for j in range (naxes[0]):
            xdist = float(j) - pos[0]
            if abs(xdist)>ignore*fwhm:
                continue
            ex = scr1[0]*xdist+scr1[1]*ydist
            ey = scr1[2]*xdist+scr1[3]*ydist
            if not dodist:
                a[i,j] = (flux/axrat)*np.exp(-(ex*ex+ey*ey))/(fwhm*fwhm*np.pi)
            else:
                a[i,j] = np.hypot(ex,ey)/1.6666667

    return a

def mxmul(a,b):
    output=np.zeros(4)
    output[0]=a[0]*b[0]+a[1]*b[2]
    output[1]=a[0]*b[1]+a[1]*b[3]
    output[2]=a[2]*b[0]+a[3]*b[2]
    output[3]=a[2]*b[1]+a[3]*b[3]
    return output

def mxinv(a):
    det=a[0]*a[3]-a[1]*a[2]
    output=np.array([a[3],-a[1],-a[2],a[0]])/det
    return output




################## closure ##############################

def closure (vis, tel, lastv=-1, plotfile='clplot.png'):
    # Find which number is which antenna
    itels = np.sort(get_idx_tels (vis, tel))
    if itels == []:
        return -1

    # Make three reference MSs with pointers
    print 'itels',itels
    d1,ut1,uvw = dget_t (vis,itels[0],itels[1])
    d2,ut2,uvw = dget_t (vis,itels[1],itels[2])
    d3,ut3,uvw = dget_t (vis,itels[0],itels[2])
    a1,p1 = getap(d1[:lastv])
    a2,p2 = getap(d2[:lastv])
    a3,p3 = getap(d3[:lastv])
    clph = p1+p2-p3
    np.putmask(clph,clph>np.pi,clph-2*np.pi)
    np.putmask(clph,clph<-np.pi,clph+2*np.pi)
    # return a statistic - is 1.64 for random closure phase, less for coherent
    if len(plotfile):
        plt.plot(clph,'b+')
        plt.savefig(plotfile)
    return np.nanmean(np.gradient(np.unwrap(clph))**2)



def closure(vis,tel,lastv=-1,pol=0,use_spw=0,bchan=0,echan=-1,doplot=False,doret=False):

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
        cp = get_amp_clph(d1[:lastv],d2[:lastv],d3[:lastv],spw[:lastv],
             pol=0, use_spw=use_spw, bchan=bchan, echan=echan)
        try:
            allcp.append(cp)
        except:
            allcp = [cp]

        os.system('rm -fr closure_temp*ms')
        if os.path.exists ('closure_which'):
            os.system('rm closure_which')
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



################### correlate ##########################

def astropy_sep (tra1,tdec1,tra2,tdec2):
    sc1 = SkyCoord(tra1, tdec1, frame = 'fk5', unit='degree')
    sc2 = SkyCoord(tra2, tdec2, frame = 'fk5', unit = 'degree')
    return((sc1.separation(sc2)).degree)

# correlate two arrays sorted in ra. array2 is the bigger one.

def correlate (array1, ra1, dec1, array2, ra2, dec2, dist, \
               mindist=0.0, isabs=False, noisy=True):
    fstart=nfstart=0
    fend=array2.shape[0]
    icou=0
    correl=np.array([])
    decfac=1.0 if isabs else min(1./np.cos(np.deg2rad(array1[:,dec1].max())),\
               1./np.cos(np.deg2rad(array2[:,dec2].max())))

    for i in range(array1.shape[0]):
        i10 = np.linspace(0,array1.shape[0],10,dtype='int')
        i100 = np.linspace(0,array1.shape[0],100,dtype='int')
        if i in i10 and noisy:
            sys.stdout.write('*')
            sys.stdout.flush()
        elif i in i100 and noisy:
            sys.stdout.write('.')
            sys.stdout.flush()
        else:
            pass
        fstart=nfstart
        for j in range(fstart,fend):
            r1,d1 = array1[i,ra1],array1[i,dec1]
            r2,d2 = array2[j,ra2],array2[j,dec2]
            radiff = r2-r1
            if radiff<-decfac*dist:
                nfstart=j
            if radiff>decfac*dist:
                break
            if abs(d2-d1)>dist:
                continue
            if isabs:
                adist = np.hypot(r1-r2,d1-d2)
            else:
		adist = astropy_sep(r1,d2,r2,d2)

            if adist<dist and abs(radiff)<90.0 and adist>=mindist:
                try:
                    correl=np.vstack((correl,np.array([i,j,adist])))
                except:
                    correl=np.array([[i,j,adist]])

    return correl



################## model_engine ########################

#--- model_engine.py: makes a (currently) two-component model
#    from u-v data on a particular closure triangle
#           v.1 Neal Jackson, 2015.09.29
#           v.2 NJ, 2017.01.16 converted to CASA, many changes
#           v.3 NJ, 2017.03.07 parallelised

## define some things
plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.interpolation']='nearest'
warnings.simplefilter('ignore')
RAD2ARC = 3600.*180./np.pi
LIGHT = 2.99792458E+8
#FIRSTNPY = './first_2008.simple.npy'
MX,MY,MF,MW,MR,MP = range(6)
ISPARALLEL = int(os.popen('nproc').read())>4

# Make a movie from a set of png files
def movie (fps=4):
    a=np.sort(glob.glob('model_engine_*.png'))
    for i in range(len(a)):
        os.system('convert %s model_engine%03d.jpeg'%(a[i],i))
    command = "mencoder \"mf://model_engine*.jpeg\" -mf fps=%d -o model_engine.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=4800" %fps
    os.system(command)
    os.system('rm model_engine*.jpeg')

# Work out real and imaginary parts of a list of visibilities given model and uvw's
def uvw2reim (uvw, model):
    u,v,re,im = uvw[:,0],uvw[:,1],0.,0.
    for m in model:     # note - sign in m[MX] as this is HA not RA
        cmpamp,cmpphs = m[MF],2.*np.pi*(-u*m[MX]+v*m[MY])/RAD2ARC
        if m[MW]!=0.0:
            sphi,cphi = np.sin(np.deg2rad(m[MP])),np.cos(np.deg2rad(m[MP]))
            tc = np.pi*(m[MW]/RAD2ARC)*np.hypot(v*cphi+u*sphi,m[MR]*(u*cphi-v*sphi))
            cmpamp = m[MF]*np.exp(-0.3696737602*tc*tc)
        re += cmpamp * np.cos(cmpphs)
        im += cmpamp * np.sin(cmpphs)
    return re,im

def dget_t (vis, tel1, tel2):
    os.system('taql \'select from %s where ANTENNA1==%d and ANTENNA2==%d giving %s\'' % (vis, tel1, tel2, 'cl_temp.ms'))
    t = pt.table('cl_temp.ms')
    ut = np.ravel(np.asarray([tuple(each.values()) for each in t.select('TIME')]))
    spw = np.ravel(np.asarray([tuple(each.values()) for each in t.select('DATA_DESC_ID')]))
    dc = t.select('DATA')
    d = np.asarray([tuple(each.values()) for each in dc])[:,0,:,:]
    d = np.swapaxes(d,0,2)     # makes the order pol - chan - time as in casa
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    if spw.sum():
        for i in np.unique(spw):
            new = np.take(d,np.argwhere(spw==i),axis=2)[:,:,:,0]
            try:
                d_out = np.concatenate((d_out,new),axis=1)
            except:
                d_out = np.copy(new)
        d = d_out
    return d,ut,uvw

def norm(a,isred=True):
    nlim = np.pi
    a = a%(2*nlim) if isred else a
    np.putmask(a,a>nlim,a-2.*nlim)
    np.putmask(a,a<-nlim,a+2.*nlim)
    return a

def getap (d,pol=0):
    ph = np.sum(d[pol],axis=0)
    return np.sum(abs(d[pol]),axis=0)/d.shape[1],np.arctan2(ph.imag,ph.real)

def get_uvw_table (t):
    for u in t.select('UVW'):
        try:
            uvw = np.vstack ((uvw,u['UVW']))
        except:
            uvw = np.copy(u['UVW'])
    return uvw

# Get data and u-v arrays from a measurement set on a given triangle
def data_extract (vis):
    global uvw01,uvw02,uvw12,cp012,a01,a02,a12
    chw = taql_calc(vis,'SPECTRAL_WINDOW','CHAN_WIDTH','mean')
    ch0 = taql_calc(vis,'SPECTRAL_WINDOW','REF_FREQUENCY','mean')
    nchan = taql_calc(vis,'SPECTRAL_WINDOW','NUM_CHAN','mean')
    sra,sdec = taql_from(vis,'FIELD','PHASE_DIR')
    nspw = taql_num(vis,'SPECTRAL_WINDOW','NUM_CHAN')
    sra = np.asarray(sra.replace('h',' ').replace('m',' ').split(),dtype='f')
    sdec = np.asarray(sdec.replace('d',' ').replace('m',' ').split(),dtype='f')
    ra = 15.0*(sra[0]+sra[1]/60.0+sra[2]/3600.0)
    dec = np.sign(sdec[0]) * ( np.abs(sdec[0])+sdec[1]/60.0+sdec[2]/3600.0 )
    wlength = np.mean(LIGHT/(ch0 + chw*np.arange(nspw*nchan)))
    itel = np.array(get_idx_tels (vis,trname))
    # order of telescopes on baselines 0-1, 0-2, 1-2; -1 if in "wrong" order
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    btel = [min(itel[0],itel[1]),max(itel[0],itel[1]),min(itel[0],itel[2]),\
            max(itel[0],itel[2]),min(itel[1],itel[2]),max(itel[1],itel[2])]
    d01,ut01,uvw01 = dget_t (vis, btel[0],btel[1])
    d02,ut02,uvw02 = dget_t (vis, btel[2],btel[3])
    d12,ut12,uvw12 = dget_t (vis, btel[4],btel[5])
    a01,p01 = getap(d01)
    a02,p02 = getap(d02)
    a12,p12 = getap(d12)
    cp012 = otel[0]*p01-otel[1]*p02+otel[2]*p12
    np.putmask(cp012,cp012>np.pi,cp012-2*np.pi)
    np.putmask(cp012,cp012<-np.pi,cp012+2*np.pi)
    uvw01 /= wlength
    uvw02 /= wlength
    uvw12 /= wlength
    print trname,'-> antenna numbers:',itel
    print 'Baseline lengths: %s-%s: %dkm %s-%s: %dkm %s-%s: %dkm' % \
       (trname[0],trname[1],int(np.sqrt((uvw01[0]**2).sum())*wlength/1000),\
        trname[0],trname[2],int(np.sqrt((uvw02[0]**2).sum())*wlength/1000),\
        trname[1],trname[2],int(np.sqrt((uvw12[0]**2).sum())*wlength/1000))
    os.system('rm -fr cl_tmp*.ms')
    return itel,np.mean(wlength),ra,dec

def model_extract (model,itel):
    otel = 1.-2.*np.asarray([itel[0]>itel[1],itel[0]>itel[2],itel[1]>itel[2]],dtype=float)
    re01,im01 = uvw2reim (uvw01,model)
    re12,im12 = uvw2reim (uvw12,model)
    re02,im02 = uvw2reim (uvw02,model)
    ph01 = norm(np.arctan2(im01,re01))
    ph02 = norm(np.arctan2(im02,re02))
    ph12 = norm(np.arctan2(im12,re12))
    clph = norm(ph01*otel[0] - ph02*otel[1] + ph12*otel[2])
    return np.hypot(re01,im01),np.hypot(re02,im02),clph

def plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou):
    ells = []
    for i in model:
        if i[MW]!=0.0:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=2.*i[MW],\
                height=2.*i[MW]*i[MR],lw=0.5,angle=i[MP]+90.,fill=0.0))
        else:
            ells.append(matplotlib.patches.Ellipse(xy=[i[MX],i[MY]],width=0.1,\
                height=0.1,lw=0.5,angle=0.,fill=1.0,color='red'))
    if plottype in [1,10]:
        plt.subplot2grid((6,5),(0,0),rowspan=2)
        scalarray = np.sort(np.ravel(a01)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,1),rowspan=2)
        plt.imshow(A01,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(0,2),rowspan=2)
        plt.imshow(a01-A01,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(2,0),rowspan=2)
        scalarray = np.sort(np.ravel(a02)); ls = len(scalarray)
        amp_min, amp_max = scalarray[int(0.1*ls)], scalarray[int(0.9*ls)]
        plt.imshow(a02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,1),rowspan=2)
        plt.imshow(A02,aspect='auto',vmin=amp_min, vmax=amp_max)
        plt.subplot2grid((6,5),(2,2),rowspan=2)
        plt.imshow(a02-A02,aspect='auto',\
                   vmin=2*amp_min-amp_max,vmax=2*amp_max-amp_min)
        plt.subplot2grid((6,5),(4,0),rowspan=2,yticks=[])
        plt.imshow(cp012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,1),rowspan=2,yticks=[])
        plt.imshow(CP012,aspect='auto',vmin=-180,vmax=180)
        plt.subplot2grid((6,5),(4,2),rowspan=2,yticks=[])
        plt.imshow(CP012-cp012,aspect='auto')
        plt.subplot2grid((6,5),(3,3),rowspan=3,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    elif plottype in [2,20]:
        plt.subplot2grid((6,5),(0,0),rowspan=2,colspan=3)
        plt.plot(a01,'b-'); plt.plot(A01,'r-')
        plt.legend([trname[0]+'-'+trname[1]],fontsize=6)
        plt.subplot2grid((6,5),(2,0),rowspan=2,colspan=3)
        plt.plot(a02,'b-'); plt.plot(A02,'r-')
        plt.legend([trname[0]+'-'+trname[2]],fontsize=6)
        plt.subplot2grid((6,5),(4,0),rowspan=2,colspan=3)
        plt.plot(cp012,'b-'); plt.plot(CP012,'r-')
        plt.subplot2grid((6,5),(4,3),rowspan=2,colspan=3,xticks=[],yticks=[])
        plt.imshow(aplot,cmap=matplotlib.cm.gray_r)
        plt.title('X2=%.3f'%goodness)
    ax = plt.subplot2grid((6,5),(0,3),rowspan=3,colspan=2)
    for i in range(len(ells)):
        e = ells[i]
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
    ax.set_xlim(-glim,glim);ax.set_ylim(-glim,glim)
    plt.grid()
    if gcou==-1:   # find how many png files, make the next one
        gcou = len(glob.glob('model_engine_*.png'))
    plt.title('model_engine_%03d.png'%gcou)
    if plottype in [10,20]:
        plt.savefig('model_engine_%03d.png'%gcou)

def ndiff (a,b):
    sqd = np.array([])
    for i in range(-len(a)/2,len(a)/2):
        a1 = np.roll(a,i)
        idx1,idx2 = max(0,i),min(len(a),len(a)+i)
        sqd = np.append(sqd, np.mean(((b-a1)[idx1:idx2]**2)))
    return sqd

def get_goodness(A01,A02,CP012):
    beta = 0.00001    #   this is a pretty vital parameter
    ascat = np.median(abs(np.gradient(np.ravel(a02))))
    cscat = np.median(abs(np.gradient(np.ravel(cp012))))
    sq = ndiff(a01*np.nanmean(A01)/np.nanmean(a01),A01)/ascat**2 + \
         ndiff(a02*np.nanmean(A02)/np.nanmean(a02),A02)/ascat**2 + \
         ndiff(cp012,CP012)/cscat**2
    difmin = 0.5*len(a01)-np.argwhere(sq==sq.min())[0][0]
    return sq.min() + beta*difmin**2

def mod_func (x0, *x):
    model,opt,itel,aplot,gcou,iy,ix = x
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = x0
    model = model.reshape(len(model)/6,6)
    A01,A02,CP012 = model_extract (model,itel)
    if ampfiddle:
        A01 *= np.median(a01)/np.median(A01)
        A02 *= np.median(a02)/np.median(A02)
    goodness = get_goodness(A01,A02,CP012)
    if plottype:
        plotimg (A01,A02,CP012,model,goodness,itel,aplot,gcou)
        if plottype in [1,2]:
            plt.draw()
            plt.pause(0.001)
            plt.clf()
    return goodness

def getmodel(coord,beam,firstnpy):
    gsiz = int(gridsize/(bsub*beam))   # original grid very big
    ginc = bsub*beam                   # arcsec/pix grid size
    grid = np.ones((gsiz,gsiz))*np.nan # original grid full of NaN
    first = np.load(firstnpy)
    pflux,pcoord = [],[]
    a = correlate (np.array([coord]),0,1,first,0,1,0.5/60)
    cosdec = np.cos(np.deg2rad(coord[1]))
    print 'Found: %d FIRST sources'%len(a)
    for i in range(len(a)):
        fi = first[int(a[i,1])]
        scoord = astropy.coordinates.SkyCoord(fi[0],fi[1],unit='degree')
        scoord = scoord.to_string(style='hmsdms')
        scoord = scoord.replace('d',':').replace('m',':')
        scoord = scoord.replace('h',':').replace('s','')
        print '%s flux=%.1f shape=(%.1fx%.1f PA%.1f)'%(scoord,fi[2],fi[4],fi[5],fi[6])
        pix = 0.5*gsiz+3600.*(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        try:
            pcoord.append(fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        except:
            pcoord = (fi[:2]-coord)*np.array([cosdec/ginc,1./ginc])
        pflux.append(fi[2])
        new = mkgauss([gsiz,gsiz],pix,fi[2],fi[4]/ginc,fi[5]/fi[4],fi[6])
        np.putmask(grid,new>0.004*new.max(),i)  # decrease 0.01 if missing cpts
    # grid now consists of original very big grid, with numbers instead of NaN where
    # the secondary might be
    if pflux:    # shrink the grid around the Gaussian near FIRST sources
        while all(np.isnan(grid[0,:])) and all(np.isnan(grid[-1,:])) and \
              all(np.isnan(grid[:,0])) and all(np.isnan(grid[:,-1])):
            grid = grid[1:-1,1:-1]
    return grid,ginc,grid.shape[0],pflux,pcoord

def parallel_function(f):
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores = max(1,int(os.popen('nproc').read())-1)  # use all cores - 1
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

def grid_search_thread (k):
    a =  mod_func(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7])
    f=open('grid_search_thread.log','a')
    f.write('%d %d %f\n'%(k[-2],k[-1],a))
    f.close()
    return a

def grid_search (model,cpt,gridcpt,itel,aplot,gcou,grid,gsiz,ginc,isparallel=ISPARALLEL):
    opt = np.zeros_like(model,dtype='bool')
    args = []
    import copy
    for ix in range(gsiz):
        x = ginc*(ix-gsiz/2.0)   # x,y in arcsec; a in ginc-size pixels
        for iy in range(gsiz):
            if grid[iy,ix] == gridcpt:
                y = ginc*(iy-gsiz/2.0)
                model[cpt][0],model[cpt][1] = x,y
                if isparallel:
                    arg = ([],model,opt,itel,aplot,gcou,iy,ix)
                    args.append(copy.deepcopy(arg))   # otherwise overwrites all elements
                    gcou += 1
                else:
                    aplot[iy][ix] = mod_func ([],model,opt,itel,aplot,gcou,iy,ix)
                    gcou += 1
    if isparallel:
        print 'Starting grid search with',len(args),'points'
        os.system('rm grid_search_thread.log')
        grid_search_thread.parallel = parallel_function(grid_search_thread)
        parallel_result = grid_search_thread.parallel (args)
        for i in range(len(parallel_result)):
            aplot[args[i][-2],args[i][-1]] = parallel_result[i]
    np.putmask(aplot,np.isnan(aplot),np.nanmax(aplot))
    model[cpt,:2] = ginc*(np.asarray(ndimage.measurements.minimum_position \
            (aplot)[::-1])-0.5*np.asarray(grid.shape))
    return model,aplot

def recentroid (model,startcpt,endcpt,ginc,gsiz):
    pos = model[startcpt:endcpt+1,:2]*model[startcpt:endcpt+1,2]
    centroid = np.sum(pos,axis=0)/len(pos)
    model[startcpt:endcpt+1,:2] -= centroid
    return model

def adjust_all (model,itel,aplot,gcou):
    opt = np.ones_like(model,dtype='bool')
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=100)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def refine_points (model,cpts,itel,aplot,gcou):
    opt = np.zeros_like(model,dtype='bool')
    for i in cpts:
        opt[i,2:] = True
    x0 = np.ravel(model)[np.ravel(opt)]
    args = (model,opt,itel,aplot,-1,0,0)
    xopt = optimize.fmin(mod_func, x0, args=args, maxiter=20)
    model,opt = np.ravel(model),np.ravel(opt)
    model[opt] = xopt
    model = model.reshape(len(model)/6,6)
    return model

def write_skymodel (ra,dec,model,outname):

    print model

    if outname!='':
        f = open(outname,'w')
        f.write ('# (Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, Orientation) = format\n')
    for i in range(len(model)):
	if isinstance( ra, str ):
	    sra = ra
	    sdec = dec
	else:
	    cosd = 3600.*np.cos(np.deg2rad(dec))
            s = astropy.coordinates.SkyCoord(ra-model[i,0]/cosd,dec+model[i,1]/3600,unit='degree')
            s = s.to_string(style='hmsdms')
            sra = s.split()[0]
	    sdec = s.split()[1]
        sra = sra.replace('h',':').replace('m',':').replace('s','')
        sdec = sdec.replace('d','.').replace('m','.').replace('s','')
        print 'ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f'%(i,sra,sdec,model[i,MF],\
              model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP]))
        if outname!='':
            f.write('ME%d, GAUSSIAN, %s, %s, %f, %f, %f, %f\n'%(i,sra,sdec,model[i,MF],\
                  model[i,MW],model[i,MW]*model[i,MR],np.rad2deg(model[i,MP])))
    if outname!='':
        f.close()

def model_engine(vis,TRNAME,firstnpy,BSUB=0.3,GRIDSIZE=12.0,PLOTTYPE=20,AMPFIDDLE=True,outname='model_engine.sky'):
    global bsub,gridsize,plottype,ampfiddle,glim,trname
    bsub,gridsize,plottype,ampfiddle,trname,gcou = BSUB,GRIDSIZE,PLOTTYPE,AMPFIDDLE,TRNAME,0
    os.system('rm model_engine*.png')
    itel,wv,ra,dec = data_extract (vis)
    s_amp = np.sort(np.ravel(a01)); ls = len(s_amp)
    flux1,flux2 = np.median(s_amp), 0.5*(s_amp[int(0.99*ls)]-s_amp[int(0.01*ls)])
    flux = np.median(a01)
    beam = RAD2ARC/np.nanmax(abs(uvw01)); print 'Beam:',beam,'arcsec'
    grid,ginc,gsiz,pflux,pcoord = getmodel (np.array([ra,dec]),beam,firstnpy)
    if not len(pflux):   # no FIRST source, search the whole grid
        grid = np.zeros_like (grid)
    aplot = np.ones_like(grid)*np.nan
    glim = 0.5*ginc*gsiz
    if plottype in [1,2]:
        plt.ion()
        fig = plt.figure()
    if len(pflux) < 2:    # zero, or one FIRST source
        model = np.array([[0.0,0.0,flux1,0.5,1.0,0.0],[-1.0,-2.0,flux2,0.5,1.0,0.0]])
        model,aplot = grid_search (model,1,0,itel,aplot,gcou,grid,gsiz,ginc)
# fix one cpt, grid-search posn of 2nd.
        print 'Model after grid search:'
        write_skymodel (ra,dec,model,'')
        model[0,MW] = model[1,MW] = 1.0
        model = refine_points (model,[0,1],itel,aplot,gcou)
# fit for size/orientation only
        print 'Model after refining points:'
        write_skymodel (ra,dec,model,'')
        model = recentroid (model,0,1,ginc,gsiz)  # put centroid in centre of image
        print 'Model after recentroiding:'
    else:   # >1 FIRST source - not tested yet
        model = np.zeros((len(pflux),6))
        for i in range(len(pflux)):
            model[i,:2] = pcoord[i]
            model[i,3] = pflux[i]*flux/pflux.sum()
            model[i,4:] = [0.5,1.0,0.0]
        model = adjust_all(model,itel,aplot,gcou)
    if plottype in [10,20]:
        movie()
    write_skymodel (ra,dec,model,'')
    write_skymodel (ra,dec,model,outname)

################## skynet ##############################

def skynet_NDPPP (vis,model,solint=1.0):
    os.system('rm -fr %s/sky\n'%vis)
    ss = 'makesourcedb in=%s out=%s/sky format=\'<\''%(model,vis)
    print ss
    os.system (ss)
    with open('NDPPP.parset','w') as f:
        f.write('msin=%s\n'%vis)
        f.write('msin.datacolumn=DATA\n')
        f.write('msout=%s\n'%vis)
        f.write('msout.datacolumn=CORRECTED_DATA\n')
        f.write('steps=[gaincal]\n')
        f.write('gaincal.applysolution=True\n')
        f.write('gaincal.sourcedb = %s/sky\n'%vis)
        f.write('gaincal.maxiter=200\n')
        f.write('gaincal.caltype=phaseonly\n')
        f.write('gaincal.solint=%d\n'%int(solint))
    f.close()
    os.system('NDPPP NDPPP.parset')


def main (vis, self_cal_script, firstnpy, mode=3, closure_tels='ST001;DE601;DE605',cthr=1.6, model_only=0, smodel=1.0, dopipe=True, lastv=-1, pol=0, use_spw=0,bchan=0,echan=-1,doplot=False,doret=False):

    ## make sure the parameters are the correct format
    mode = int( mode )
    cthr = float( cthr )
    smodel = float( smodel )
    model_only = int( model_only )

    if not len(closure_tels) == 3:
	closure_tels = closure_tels.split(';')
    
    ra,dec = taql_from (vis, 'FIELD', 'PHASE_DIR')
    
    ## prepare the information for getting the closure phase scatter
    tel = np.array(closure_tels.split(';'))
    if len(tel)%3:
        print 'Station list is not a multiple of 3 long, truncating it...'
        tel = tel[0:3*(len(tel)//3)]
    if dopipe and len(tel)>3:
        print 'dopipe True, truncating to 3 telescopes'
        tel = tel[:3]
    tel = tel.reshape(len(tel)/3,3)
    closure_scatter = closure(vis, tel,lastv=lastv,pol=pol,use_spw=use_spw,bchan=bchan,echan=echan,doplot=doplot,doret=doret)
    if dopipe:
	print '\n Scatter for the direction ' + vis.split('/')[-1].split('_')[0] + ' is %s \n' % closure_scatter
        os.system('echo Scatter for the direction ' + vis.split('/')[-1].split('_')[0] + ' is ' + str(scatter_cp) + ' >> closure_phases.txt')
    else:
	try:
	    allret.append(closure_scatter)
	except:
	    allret = [closure_scatter]

    ## check if it's sensible
    if dopipe:
        if closure_scatter > cthr:
            return closure_scatter
    if mode == 1:     # this is the default of the self_calibration_pipeline_V2.
	print 'mode 1: wsclean model'
        os.system('python '+self_cal_script+' '+vis+' -p') # amplitudes?
    if mode == 2:   # use NDPPP to selfcal the long baselines to a small (0.1") Gaussian
	print 'mode 2: point model'
	point_model = np.array( [ [0.0,0.0,smodel,0.1,0.0,0.0] ] )
	write_skymodel (ra,dec,point_model,vis+'/skymodel')
	if model_only == 0:
            skynet_NDPPP (vis,vis+'_mod',solint=5)  # timesteps
            os.system('python '+self_cal_script+' -d CORRECTED_DATA '+vis+' -p')
	else:
	    # run makesourcedb to generate sky
	    print 'MODEL ONLY'
	    ss = 'rm -fr %s/sky\n'%vis
	    print ss
	    os.system( ss )
	    ss = 'makesourcedb in=%s/skymodel out=%s/sky format=\'<\''%(vis,vis)
	    print ss
	    os.system ( ss )
    if mode == 3:   # make an engine model and selfcal against this
	print 'mode 3: model_engine model'
        model_engine (vis,closure_tels,firstnpy,PLOTTYPE=0,outname=vis+'_mod')
	if model_only == 0:
            skynet_NDPPP (vis,vis+'_mod',solint=5)
            os.system('python '+self_cal_script+' -d CORRECTED_DATA '+vis+' -p')
    #else:
    #    return np.nan
    #return 0.0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Skynet script to handle LBCS calibrators.")

    parser.add_argument('vis', type=str, help='Measurement set for which to run skynet')
    parser.add_argument('--self_cal_script',type=str, help='Self-calibration script to use')
    parser.add_argument('--firstnpy',type=str,help='absolute path of first_2008.simple.npy')
    parser.add_argument('--mode',type=int, help='Mode to use', default=3)
    parser.add_argument('--closure_tels',type=str,help='Stations to use for calculating closure phase.', default='ST001;DE601;DE605' )
    parser.add_argument('--cthr',type=float,help='Threshold for closure phase scatter.', default=1.6)
    parser.add_argument('--smodel',type=float,help='Flux density in Jy of point source.', default=1.0)
    parser.add_argument('--model_only',type=int,help='set to 1 to get model only',default=0)

    args = parser.parse_args()

    main( args.vis, args.self_cal_script, args.firstnpy, mode=args.mode, closure_tels=args.closure_tels, cthr=args.cthr, smodel=args.smodel, model_only=args.model_only )

