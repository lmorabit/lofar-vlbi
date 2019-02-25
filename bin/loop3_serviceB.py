# loop3_service   Neal Jackson 1.10.2018
# Service routines for loop3. 
# Response to complaints that this looks too much like AIPS will involve adding 
#    APARM arrays as arguments.

import numpy as np,os,sys,glob,time,scipy,pickle
import losoto.h5parm as h5parm
from scipy import interpolate
import pyrap.tables as pt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import h5py

def loop3log (vis, pstr, cret = True):
    fo = open(vis+'_proc.log','a')
    fo.write('%s'%pstr)
    fo.write('\n' if cret else '')
    fo.close()
    print pstr

# h5 routine to read h5 files. This is done in a separate python call
# otherwise the hparm does not close properly. No info found on any other
# way to do this.

def h5read (htab, solset, soltab):
    fo = open('tmp_read.py','w')
    fo.write ('import losoto.h5parm as h5parm,os,pickle\n')
    fo.write ('tab = h5parm.openSoltab(\'%s\',solsetName=\'%s\',soltabName=\'%s\')\n' % \
                          (htab,solset,soltab)   )
    fo.write ('v, vm = tab.getValues()[0], tab.getValues()[1]\n')
    fo.write ('pickle.dump(v,open(\'v.pkl\',\'wb\'))\n')
    fo.write ('pickle.dump(vm,open(\'vm.pkl\',\'wb\'))\n')
    fo.close()
    os.system ('python tmp_read.py')
    v = pickle.load (open('v.pkl','rb'))
    vm = pickle.load (open('vm.pkl','rb'))
    os.system('rm tmp_read.py;rm v.pkl;rm vm.pkl')
    return v, vm

# Return a solution to zero (and ones if an amplitude solution exists).
# Used after detection of an incoherent solution on an antenna.
def zerosol (H1,ant):
    h1 = h5py.File(H1,'r+')
    n1 = h1.get('sol000/phase000')
    v1 = np.array(n1['val'])
    z = np.zeros_like(v1[:,:,0,:])
    ant1 = np.array(h1.get('sol000/phase000/ant'))
    for i in range(len(ant1)):
        if ant1[i] in ant:
            try:
                h1['sol000/amplitude000/val'][:,:,i,:] = z+1.0
            except:
                pass
            h1['sol000/phase000/val'][:,:,i,:] = z
    h1.close()

# Given two calibration structures H1 and H2, and antennas to interpolate, 
# replace the phase calibration of H1 for each antenna with an interpolated 
# version of H2. (Has been tested for phase, needs testing for amplitude)
def clcal (H1,H2,ant_interp=None):
    isamp = True
    h1,h2 = h5py.File(H1,'r+'),h5py.File(H2)
    n1,n2 = h1.get('sol000/phase000'),h2.get('sol000/phase000')
    t1,t2 = np.array(n1['time']),np.array(n2['time'])
    v1,v2 = np.array(n1['val']),np.array(n2['val'])
    a1 = np.array(h1.get('sol000/phase000/ant'))
    try:
        na1,na2 = h1.get('sol000/amplitude000'),h2.get('sol000/amplitude000')
        va1,va2 = np.array(na1['val']),np.array(na2['val'])
    except:
        isamp = False
    ant_interp = a1 if ant_interp==None else ant_interp
    for i in range(len(a1)):
        # SM: Crashed here for 3C 273 (but not for 3C 280). Traceback:
        # File "/data020/scratch/sean/fafa/lofar-lb/loop3_serviceB.py", line 72, in clcal
        # for i in range(len(a1)):
        # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()

        if a1[i] not in ant_interp:
            continue
        for iz in range(v1.shape[1]):
            for ipol in range(v1.shape[3]):
                if isamp:
                    zr = va2[:,iz,i,ipol]*np.cos(v2[:,iz,i,ipol])
                    zc = va2[:,iz,i,ipol]*np.sin(v2[:,iz,i,ipol])
                    z2 = zr + 1j*zc
                    z = scipy.interpolate.griddata(t2,z2,t1,method='linear')
                    h1['sol000/amplitude000/val'][:,iz,i,ipol] = abs(z)
                    h1['sol000/phase/val'][:,iz,i,ipol] = np.arctan2(z.imag,z.real)
                else:
                    z = scipy.interpolate.griddata(t2,\
                        np.unwrap(v2[:,iz,i,ipol]),t1,method='linear')
                    while z.max()>np.pi:
                        np.putmask(z,z>np.pi,z-2.*np.pi)
                    while z.min()<-np.pi:
                        np.putmask(z,z<-np.pi,z+2.*np.pi)
                    h1['sol000/phase000/val'][:,iz,i,ipol] = z

    h1.close(); h2.close()
    h1 = h5py.File(H1,'r+')
    n1 = h1.get('sol000/phase000')
    v1 = np.array(n1['val'])
    h1.close()

def calib (vis,incol='DATA',outcol='DATA',solint=180,solmode='P',\
           model=None,outms='.',outcal=None,tsamp=8.0):
    loop3log (vis,'-------> %d %.1f'%(solint,tsamp))
    outcal = vis+'_cal' if outcal==None else outcal
    mgain = 'sourcedb=%s\n'%model if model else 'usemodelcolumn=true\n'
    caltype = 'phaseonly' if solmode=='P' else 'diagonal'
    f=open('calib.parset','w')
    f.write('msin=%s\n'%vis)
    f.write('msin.datacolumn=%s\n'%incol)
    f.write('msout=%s\n'%outms)
    f.write('msout.datacolumn=%s\n'%outcol)
    f.write('steps=[gaincal]\n')
    f.write('gaincal.'+mgain)
    f.write('gaincal.caltype=%s\n'%caltype)
    f.write('gaincal.solint=%i\n'%(solint/tsamp))
    f.write('gaincal.usebeammodel=False\n')
    f.write('gaincal.parmdb=%s\n'%outcal)
    f.write('gaincal.applysolution=%s\n'%('False' if incol==outcol else 'True'))
    f.close()
    time_start = time.time()
    # Bug fix here: NDPPP leaves the .h5 files unclosed. So we have to 
    # start a separate python session to run the NDPPP on calib.parset, 
    # which closes the .h5 files on exit.
    fo=open('calib.py','w')
    fo.write ('import os\nos.system(\'NDPPP calib.parset\')\n')
    fo.close()
    os.system('python calib.py')
    time_end = time.time()
    loop3log(vis,'NDPPP took %d s' % int(time_end-time_start))


# Make the coherence parameter. This relies on the difference in the phase
# solutions in XX and YY remaining constant if the solutions are coherent.
# Also need to return an incoherent answer (2.0) if there are too many NaN
# solutions (here >10%)
def coherence_metric (htab='1327_test.ms_cal.h5',solset='sol000',soltab='phase000'):
    NANFRAC, INCOH = 0.1, 2.0
    v, vm = h5read (htab, solset, soltab)
    ant,freq,pol,time = vm['ant'],vm['freq'],vm['pol'],vm['time']
    coh = np.array([])
    for i in range(len(ant)):   # assumes two polarizations XX YY
#     changed this (njj) - note that np.unwrap gives an array full of NaN
#     if even the first element of the input array is NaN
#        diff = np.unwrap(v[:,0,i,0]-v[:,0,i,1])
        diff = v[:,0,i,0]-v[:,0,i,1]
        if float(len(diff[np.isnan(diff)]))>NANFRAC*float(len(diff)):
            coh = np.append(coh,INCOH)
        else:
            diff = np.unwrap(diff[~np.isnan(diff)])
            coh = np.append(coh,np.nanmean(np.gradient(abs(diff))**2))
    return coh
    

def snplt (vis,htab='1327_test.ms_cal.h5',solset='sol000',soltab='phase000',\
           antenna=None,nplot=6,outpng=None):
    outpng = outpng if outpng else htab
    v,vm = h5read (htab, solset, soltab)
    ant,freq,pol,time = vm['ant'],vm['freq'],vm['pol'],vm['time']
    time = 24.*(time/86400. - int(time[0])/86400)
    iplot = 0
    antenna = antenna if antenna else ant
    plt.clf()
    while iplot<len(antenna):
        a = antenna[iplot]
        aidx = np.argwhere(ant==a)[0][0]
        sys.stdout.write(a+' ')
        for ipol in range(v.shape[3]):
            if not (iplot+1)%nplot:
                plt.subplot(nplot,1,1+iplot%nplot)
            else:
                plt.subplot(nplot,1,1+iplot%nplot,xticks=[])
            if soltab[:5]=='phase':
                plt.plot(time,np.rad2deg(v[:,0,aidx,ipol]),'+')
                try:
                    plt.ylim(-180.,180.);plt.xlim(time[0],time[-1])
                except:
                    print 'Could not set one of the x axis limits. One of these is not real:', time[0], time[-1]
                    print 'Hard-coding to plt.xlim(0, 600)'
                    plt.ylim(-180.,180.);plt.xlim(0, 600)
                plt.text(time[0],180.0-12.*nplot,a)
            else:
                plt.plot(time,v[:,0,aidx,ipol],'+')
                vmin,vmax = min(v[:,0,aidx,ipol]),max(v[:,0,aidx,ipol])
                try:
                    plt.ylim(vmin,vmax);plt.xlim(time[0],time[-1])
                except:
                    print 'Could not set one of the axes limits. One of these is not real:', vmin, vmax, time[0], time[-1]
                plt.text(time[0],vmin+0.9*(vmax-vmin),a)
            plt.subplots_adjust(wspace=0,hspace=0)
        iplot+=1
        if not iplot%nplot:
            thispng = outpng+'_%d.png'%(iplot//nplot -1)
            if os.path.isfile(thispng):
                os.system('rm %s'%thispng)
            loop3log(vis,'-> %s'%thispng)
            try:
                plt.savefig(thispng) # ,bbox_inches='tight')
            except:
                print 'Failed to save', thispng
            plt.clf()
    if iplot%nplot:
        thispng = outpng+'_%d.png'%(iplot//nplot)
        if os.path.isfile(thispng):
            os.system('rm %s'%thispng)
        loop3log(vis,'-> %s'%thispng)
        try:
            plt.savefig(thispng) #,bbox_inches='tight')
        except:
            print 'Failed to save', thispng


# Because I don't like writing enormous command lines in code. Also only have to change once if the
# wsclean arguments change - or indeed if we use a different imager.
def imagr (vis,threads=0,mem=100,doupdatemodel=True,tempdir='',dosaveweights=False,doprimary=False,\
           robust=-1,domfsweight=False,gausstaper=0.0,tukeytaper=0.0,dostoreweights=False,outname='wsclean',\
           imsize=1024,cellsize='0.05asec',dopredict=False,niter=10000,pol='I',datacolumn='',autothreshold=3.,\
           dolocalrms=False,gain=0.1,mgain=1.0,domultiscale=False,dojoinchannels=False,channelsout=0,fitsmask='',\
           baselineaveraging=0.0,maxuvwm=0.0,minuvwm=100000.0,maxuvl=0.0,minuvl=0.0,dostopnegative=False,automask=0.,\
           dosavesourcelist=False,weightingrankfilter=0.0,weightingrankfiltersize=0.0):
    cmd = 'wsclean '
    cmd += ('' if not threads else '-j '+str(threads)+' ')
    cmd += ('' if not mem==100 else '-mem '+str(mem)+' ')
    cmd += ('' if doupdatemodel else '-no-update-model-required ')
    cmd += tempdir+' '
    cmd += ('' if not dosaveweights else '-save-weights ')
    cmd += ('' if not doprimary else '-apply-primary-beam ')
    if robust >=5:
        cmd += '-weight natural '
    elif robust <=-5:
        cmd += '-weight uniform '
    else:
        cmd += '-weight briggs %f '%robust
    cmd += ('' if not domfsweight else '-mfs-weighting ')
    cmd += ('' if gausstaper==0.0 else '-taper-gaussian %f '%gausstaper)
    cmd += ('' if tukeytaper==0.0 else '-taper-tukey %f '%tukeytaper)
    cmd += ('' if not dostoreweights else '-store-imaging-weights ')
    cmd += '-name '+outname+' '
    cmd += '-size '+str(imsize)+' '+str(imsize)+' '
    cmd += '-scale '+str(cellsize)+' '
    cmd += ('' if not dopredict else '-predict ')
    cmd += ('-niter '+str(niter)+' ')
    cmd += ('' if pol=='I' else '-pol '+pol+' ')
    cmd += ('' if datacolumn=='' else '-datacolumn %s '%datacolumn)
    cmd += ('' if autothreshold==0. else '-auto-threshold %f '%autothreshold)
    cmd += ('' if not dolocalrms else '-local-rms ')
    cmd += ('' if not domultiscale else '-multiscale ')
    cmd += ('' if channelsout==0 else '-channels-out %d -join-channels '%channelsout )    
    cmd += ('' if gain==0.1 else '-gain %f '%gain)
    cmd += ('' if mgain==1.0 else '-mgain %f '%mgain)
    cmd += ('' if fitsmask=='' else '-fits-mask %s '%fitsmask)
    cmd += ('' if baselineaveraging==0.0 else '-baseline-averaging %f '%baselineaveraging)
    cmd += ('' if maxuvl==0.0 else '-maxuv-l %f '%maxuvl)
    cmd += ('' if minuvl==0.0 else '-minuv-l %f '%minuvl)
    cmd += ('' if maxuvwm==0.0 else '-maxuvw-m %f '%maxuvwm)
    cmd += ('' if minuvwm==0.0 else '-minuvw-m %f '%minuvwm)
    cmd += ('' if not dostopnegative else '-stop-negative ')
    cmd += ('' if automask==0. else '-auto-mask %f '%automask)
    cmd += ('' if not dosavesourcelist else '-save-source-list ')
    cmd += ('' if weightingrankfilter==0.0 else '-weighting-rank-filter %f '%weightingrankfilter)
    cmd += ('' if weightingrankfiltersize==0.0 else '-weighting-rank-filter-size %f '%weightingrankfiltersize)
    cmd += vis+ '>>wsclean_chunterings'
    loop3log (vis,'Executing: '+cmd)
    os.system (cmd)


def getcoh_baseline (antenna_list, coh, ccrit):
    '''Returns the maximum baseline length for imaging given
    a list of coherences for stations'''
    aname = ['DE601','DE602','DE603','DE604','DE605','DE609','SE','FR',\
             'UK','PL','IE']
    alen = [260,580,400,420,230,200,600,700,602,800,800]
    cohlength = 2000.0
    np.putmask(coh,coh==-1.0,ccrit)
    for i in range(len(antenna_list)):
        for j in range(len(aname)):
            if antenna_list[i][:len(aname[j])]==aname[j]:
                if coh[i]>ccrit-0.1:
                    cohlength = alen[j]
                    break
    return 1000.0*cohlength
        

def montage_plot(vis):
    import glob
    imgroot = np.sort(glob.glob(vis+'*MFS-image.fits'))
    nloop = len(imgroot)
    h5png = glob.glob(vis+'*h5_*.png')
    npng = 0
    for i in h5png:
        npng = max(npng,int(i.split('_')[-1].split('.')[0])+1)
    cmd = 'montage -tile %dx%d -geometry 600x600 '%(npng+1,nloop)
    for i in range(nloop):
        os.system('python /data020/scratch/sean/fafa/lofar-lb/aplpy_makeplot.py '+imgroot[i])  # otherwise it does not find the script
        thisv = imgroot[i].split('-MFS-image.fits')[0]
        for j in range(npng):
            this = '%s_c0.h5_%d.png'%(thisv,j)
            cmd += (this+' ') if os.path.isfile(this) else 'null: '
        this = thisv+'-MFS-image.png'
#        print 'trying to add',this,os.path.isfile(this)
        cmd += (this+' ') if os.path.isfile(this) else 'null: '
    cmd += '%s_output.png'%vis
    print cmd
    os.system(cmd)
