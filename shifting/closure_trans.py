import Wizardry, numpy as np
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from AIPS import AIPS, AIPSDisk; from AIPSTask import AIPSTask, AIPSList
import matplotlib; from matplotlib import pyplot as plt

def closure_aips (aipsno, triangle, inna, incl, indisk=1, inseq=1, minmatch=5, ipol=0, dt=0.5, spw=0, plotfile=''):
    AIPS.userno = aipsno
    data = WizAIPSUVData (inna, incl, indisk, inseq)
    antnum = np.array([],dtype='int')
    for i in range(3):
        foundit = False
        for j in range (len(data.antennas)):
            tellen = min(len(data.antennas[j]),len(triangle[i]))
            if tellen<minmatch:
                if triangle[i][:tellen] == data.antennas[j][:tellen]:
                    foundit = True
                    antnum = np.append(antnum,j+1)
                    break
            else:
                if triangle[i][:minmatch] == data.antennas[j][:minmatch]:
                    foundit = True
                    antnum = np.append(antnum,j+1)
                    break
                    
    antnum = np.sort(antnum)
    isbase = [[antnum[0],antnum[1]],[antnum[1],antnum[2]],[antnum[0],antnum[2]]]
    icou, nvis = 0, len(data)
    d = [np.array([]),np.array([]),np.array([])]
    t = [np.array([]),np.array([]),np.array([])]
    for v in data:
        if v.baseline in isbase:
            vwhich = isbase.index(v.baseline)
        else:
            continue        
        vvis = v.visibility
        vv = vvis.reshape(vvis.shape[0]*vvis.shape[1],vvis.shape[2],vvis.shape[3])
        vphas = np.arctan2(vv[:,ipol,1],vv[:,ipol,0])
        np.putmask(vphas,vv[:,ipol,2]==0.0,np.nan)
        try:
            d[vwhich] = np.vstack((d[vwhich],vphas))
            t[vwhich] = np.append(t[vwhich],86400.*v.time)
        except:
            d[vwhich] = np.copy(vphas)
            t[vwhich] = np.append(t[vwhich],86400.*v.time)
        icou+=1
        if icou%100000==0:
            print 'Visibility %d/%d' % (icou,nvis)

    tu = np.sort(np.unique(np.append(t[0],np.append(t[1],t[2]))))
    for i in range(2):
        tu_diff = np.append(np.array([1.0E9]),np.diff(tu))
        tu = tu[tu_diff>dt]

#    print len(t[0]),len(t[1]),len(t[2]),'********'
    for i in range(3):
        if len(t[i])==len(tu):
            continue
        j=0; k=0; kmiss=[]
        dummy = np.ones_like(d[i][0])*np.nan
        while j<len(tu) and k<len(t[i]):
            if abs(tu[j]-t[i][k])<dt:
                j+=1
                k+=1
            else:
                kmiss.append(k)
                d[i] = np.insert (d[i],j+1,dummy,axis=0)
                j+=1
        while len(d[i])<len(tu):       # missing data at the end
            d[i] = np.insert(d[i],len(d[i]),dummy,axis=0)

    np.save('d0_aips',d[0]);np.save('d1_aips',d[1]);np.save('d2_aips',d[2])
    from scipy.fftpack import *
    p = np.asarray(d[0]+d[1]-d[2],dtype='float')
    np.putmask(p,p>np.pi,p-2.*np.pi)   # no need to worry about nans
    np.putmask(p,p<-np.pi,p+2.*np.pi)
    prand = np.random.random(p.shape[0]*p.shape[1]).reshape(p.shape)*2*np.pi-np.pi
    np.putmask(p,np.isnan(p),prand)
    pp = np.cos(p)+1j*np.sin(p)
    pprand = np.cos(prand)+1j*np.sin(prand)
    qbig = fftshift(fft2(pp))
    qbigrand = fftshift(fft2(pprand))
    yr = (qbig.shape[0]/2 - 30,qbig.shape[0]/2+30)
    xr = (qbig.shape[1]/2 - 30,qbig.shape[1]/2+30)
    q = qbig[yr[0]:yr[1],xr[0]:xr[1]]
    qrand = qbigrand[yr[0]:yr[1],xr[0]:xr[1]]
    qq = np.sort(np.ravel(abs(q)))
    qqrand = np.sort(np.ravel(abs(qrand)))
    x,y = np.arange(-30,30),np.arange(-30,30)
    xx,yy = np.meshgrid(x,y)
    z = np.hypot(xx,yy)
    maxsig = 0.0
    for i in range(3,30):
        zz=np.zeros_like(z)
        np.putmask(zz,z<i,1.0)
        sig = ((zz*abs(q)).sum()-(zz*abs(qrand)).sum()) \
                  /(np.std(abs(qrand))*np.sqrt(zz.sum()))
        maxsig = max(maxsig,sig)
    if plotfile!='':
        plt.subplot(111,xticks=[],yticks=[])
        plt.imshow(abs(q))
        plt.title(plotfile)
        plt.text(5,5,'%.1f'%maxsig)
        plt.savefig(plotfile,bbox_inches='tight')
        plt.clf()
    return maxsig

def phavg (phase, n):    # average an array of phases in chunks of size n
    phase = phase[0:n*(len(phase)//n)]
    rreal,rimag = np.cos(phase), np.sin(phase)
    arreal = np.average(np.reshape(rreal,(-1,n)),axis=1)
    arimag = np.average(np.reshape(rimag,(-1,n)),axis=1)
    return np.arctan2(arimag,arreal)
